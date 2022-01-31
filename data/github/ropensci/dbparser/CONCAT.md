
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
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit attributes file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
library(dplyr)
library(ggplot2)
```
# dbparser <img src="man/figures/logo.png" align="right" />

[![Build Status](https://travis-ci.org/ropensci/dbparser.svg?branch=master)](https://travis-ci.org/ropensci/dbparser) 
[![Build status](https://ci.appveyor.com/api/projects/status/k18sqp55n39f3y5w?svg=true)](https://ci.appveyor.com/project/MohammedFCIS/dbparser)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/dbparser)](https://cran.r-project.org/package=dbparser)
[![codecov](https://codecov.io/gh/ropensci/dbparser/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/dbparser)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://www.tidyverse.org/lifecycle/#stable)
[![metacran downloads](https://cranlogs.r-pkg.org/badges/grand-total/dbparser)](https://cran.r-project.org/package=dbparser)
[![Rdoc](https://www.rdocumentation.org/badges/version/dbparser)](https://www.rdocumentation.org/packages/dbparser)
[![CII Best Practices](https://bestpractices.coreinfrastructure.org/projects/3311/badge)](https://bestpractices.coreinfrastructure.org/projects/3311)
[![](https://badges.ropensci.org/347_status.svg)](https://github.com/ropensci/software-review/issues/347)

## Introduction

[DrugBank](https://www.drugbank.ca/) database is a comprehensive, freely
accessible, online database containing information on drugs and drug
targets. As both a bioinformatics and a cheminformatics resource,
DrugBank combines detailed drug (i.e. chemical, pharmacological and
pharmaceutical) data with comprehensive drug target (i.e. sequence,
structure, and pathway) information. More information about DrugBank can
be found [here](https://www.drugbank.ca/about).

In its raw form, the DrugBank database is a single
XML file. Users must create an [account](https://www.drugbank.ca/public_users/sign_up)
with DrugBank and request permission to [download](https://www.drugbank.ca/releases/latest)
the database. Note that this may take a couple of days.

The `dbparser` package parses the DrugBank XML database into `R` tibbles that can be explored and analyzed by the user, check [this tutorial](https://docs.ropensci.org/dbparser/articles/dbparser.html) for more details.

Also, the package offers the option to save these tibbles in databases including **SQL Server DB** and **Maria DB** just by enabling `save_table` option,  check [this tutorial](https://docs.ropensci.org/dbparser/articles/Database_Saving.html) for more details.

If you are waiting for access to the DrugBank database, or do not intend to do a deep dive with
the data, you may wish to use the `dbdataset`
[package](https://mohammedfcis.github.io/dbdataset/index.html), which contains
the DrugBank database already parsed into `R` tibbles. Note that this is a large package that
exceeds the limit set by CRAN. It is only available on GitHub.

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
Please note that the 'dbparser' project is released with a
  [Contributor Code of Conduct](https://docs.ropensci.org/dbparser/CODE_OF_CONDUCT.html).
  By contributing to this project, you agree to abide by its terms.
  
## Contributing Guide
👍🎉 First off, thanks for taking the time to contribute! 🎉👍
Please review our [Contributing Guide](https://docs.ropensci.org/dbparser/CONTRIBUTING.html).

## Share the love ❤️

Think **dbparser** is useful? Let others discover it, by telling them in person, via Twitter or a blog post.

Using **dbparser** for a paper you are writing? Consider citing it
```{r}
citation("dbparser")
```
---
title: "DrugBank Database XML Parser"
author: "Mohammed Ali, Ali Ezzat"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{DrugBank Parser}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "docs/articles/",
  out.width = "100%"
)
```

## Introduction
The main purpose of the `dbparser` package is to parse the 
[DrugBank](https://www.drugbank.ca/) database which is downloadable in XML format 
from [this link](https://www.drugbank.ca/releases/latest). The parsed data can 
then be explored and analyzed as desired by the user. 
In this tutorial, we will see how to use `dbparser` along with `dplyr` and 
`ggplot2` along with other libraries to do simple drug analysis

## Loading and Parsing the Data

Before starting the code we are assuming the following:

- user already downloaded *DrugBank* xml database file based on the
[Read Me](https://docs.ropensci.org/dbparser/) instructions or the above note,
- user saved the downloaded database in working directory as `C:\`.
- user named the downloaded xml file **drugbank.xml**. 

Now we can loads the `drugs` info, `drug groups` info and `drug targets`
actions info.

```{r eval=FALSE}
## load dbparser package
library(dbparser)
library(dplyr)
library(ggplot2)
library(XML)

## parse data from XML and save it to memory
read_drugbank_xml_db("C:\drugbank.xml")

## load drugs data
drugs <- drug()

## load drug groups data
drug_groups <- drug_groups()

## load drug targets actions data
drug_targets_actions <- targets_actions()
```
### Saving into Database
User can save parsed tibbles in database by enabling `save_table` option. More details about that can be found in [database tutorial](https://docs.ropensci.org/dbparser/articles/Database_Saving.html)

## Exploring the data

Following is an example involving a quick look at a few aspects of the parsed 
data. First we look at the proportions of `biotech` and `small-molecule` drugs 
in the data.

```{r eval=FALSE}
## view proportions of the different drug types (biotech vs. small molecule)
drugs %>% 
    select(type) %>% 
    ggplot(aes(x = type, fill = type)) + 
    geom_bar() + 
    guides(fill = FALSE)     ## removes legend for the bar colors
```


Below, we view the different `drug_groups` in the data and how prevalent they 
are.

```{r eval=FALSE}
## view proportions of the different drug types for each drug group
drugs %>% 
    full_join(drug_groups, by = c('primary_key' = 'drugbank_id')) %>% 
    select(type, group) %>% 
    ggplot(aes(x = group, fill = type)) + 
    geom_bar() + 
    theme(legend.position = 'bottom') + 
    labs(x = 'Drug Group', 
         y = 'Quantity', 
         title = "Drug Type Distribution per Drug Group", 
         caption = "created by ggplot") + 
    coord_flip()
```

Finally, we look at the `drug_targets_actions` to observe their proportions as 
well.

![](fig2.png)

```{r eval=FALSE}
## get counts of the different target actions in the data
targetActionCounts <- 
    drug_targets_actions %>% 
    group_by(action) %>% 
    summarise(count = n()) %>% 
    arrange(desc(count))

## get bar chart of the 10 most occurring target actions in the data
p <- 
    ggplot(targetActionCounts[1:10,], 
           aes(x = reorder(action,count), y = count, fill = letters[1:10])) + 
    geom_bar(stat = 'identity') +
    labs(fill = 'action', 
         x = 'Target Action', 
         y = 'Quantity', 
         title = 'Target Actions Distribution', 
         subtitle = 'Distribution of Target Actions in the Data',
         caption = 'created by ggplot') + 
    guides(fill = FALSE) +    ## removes legend for the bar colors
    coord_flip()              ## switches the X and Y axes

## display plot
p
```

![](fig3.png)
---
title: "Saving data into Database"
author: "Mohammed Ali"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Saving Data into Database}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "docs/articles/",
  out.width = "100%"
)
```

## Introduction
This tutorial aims to explain how `dbparser` can work along with R database functionalities to save parsed drug tibbles to the user desired databases.
This tutorial addresses the following three options:

- [SQLite](https://github.com/r-dbi/RSQLite)
- [RDBMS](https://db.rstudio.com/dbi/)
- [Maria Knowledge Base](https://github.com/r-dbi/RMariaDB)


*Please note that this tutorial does not explain how to install these databases as it is out of scope.*

### SQLite
SQLite is an inmemory database you can use locally easily. To save drug information using this database run the following
```{r eval = FALSE}
# Load dbparser package
library(dbparser)
# Create SQLite database connection
database_connection <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
# DrugBank database sample name
biotech <- "drugbank_record_biotech.xml"
# Use DrugBank sample database in the library
read_drugbank_xml_db(system.file("extdata", biotech, package = "dbparser"))
# Parse all available drug tibbles
run_all_parsers(save_table = TRUE, database_connection = database_connection)
# List saved tables
DBI::dbListTables(database_connection)
# Close SQLite connection
DBI::dbDisconnect(database_connection)
```


### DBI Supported Databases

`DBI` separates the connectivity to the DBMS into a *front-end* and a *back-end*. Applications use only the exposed front-end API. The back-end facilities that communicate with specific DBMSs (SQLite, MySQL, PostgreSQL, MonetDB, etc.) are provided by drivers (other packages) that get invoked automatically through S4 methods.
For more information about DBI package please refer to [this link](https://db.rstudio.com/dbi/)

The following are two examples of how to make the connection with *SQL Server* and *Maria DB*

#### SQL Server

* Make sure you have a working connection to *SQL Server* instance
![](sql_connection.png)
* Create new empty database to store drug information
![](new_database.png)
* Execute the following commands
```{r eval=FALSE}
# Load dneeded packages
library(dbparser)
library(odbc)
# Create SQLServer database connection
con <- DBI::dbConnect(odbc::odbc(), Driver = "SQL Server", Server = "MOHAMMED\\SQL2016", 
    Database = "drugbank", Trusted_Connection = T)
# Use DrugBank sample database in the library
biotech <- "drugbank_record_biotech.xml"	
# DrugBank database sample name
read_drugbank_xml_db(system.file("extdata", biotech, package = "dbparser"))
# Parse all available drug tibbles
run_all_parsers(save_table = TRUE, database_connection = con)
# List saved tables
DBI::dbListTables(con)
# Close SQLServer connection
DBI::dbDisconnect(con)
```

Then refresh your database to see new tables
![](sql_tables.png)

#### Maria Knowledge Base

* Install MariaDB
* Install MySQL client for MariaDB
* Create new empty database to store drug information

![](new_database_2.png)

* Execute the following commands
```{r eval=FALSE}
# Load dneeded packages
library(dbparser)
library(RMariaDB)
# Create SQLServer database connection
con <- RMariaDB::dbConnect(RMariaDB::MariaDB(), Server = "MariaDB",
                           dbname = "drugbank", username="root",
                           password="root")
# Use DrugBank sample database in the library
biotech <- "drugbank_record_biotech.xml"	
# DrugBank database sample name
read_drugbank_xml_db(system.file("extdata", biotech, package = "dbparser"))
# Parse all available drug tibbles
run_all_parsers(save_table = TRUE, database_connection = con)
# List saved tables
RMariaDB::dbListTables(con)
# Close SQLServer connection
RMariaDB::dbDisconnect(con)
```

Then refresh your database to see new tables
![](maria_tables.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_categories}
\alias{drug_categories}
\title{Drug Categories parser}
\usage{
drug_categories(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 2 variables:
\describe{
 \item{category}{category name}
 \item{mesh-id}{The Medical Subjects Headings (MeSH) identifier for the
 category.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
General categorizations of the drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_ahfs_codes}
\alias{drug_ahfs_codes}
\title{Drug ahfs-codes parser}
\usage{
drug_ahfs_codes(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{ahfs-code}{}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
The American Hospital Formulary Service (AHFS) identifier for this drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_food_interactions}
\alias{drug_food_interactions}
\title{Drug Groups parser}
\usage{
drug_food_interactions(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{food-interaction}{}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Food that may interact with this drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_mixtures}
\alias{drug_mixtures}
\title{Drug Mixtures parser}
\usage{
drug_mixtures(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 4 variables:
\describe{
 \item{name}{The proprietary name provided by the manufacturer for this
 combination product.}
 \item{ingredients}{A list of ingredients, separated by addition symbols}
 \item{supplemental-ingredients}{List of additional active ingredients which
  are not clinically relevant to the main indication of the product,
  separated by addition symbols.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
All commercially available products in which this drug is available in
combination with other drug molecules
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_pathway_node_parsers.R
\name{drug_pathway}
\alias{drug_pathway}
\title{Drug Pathway parser}
\usage{
drug_pathway(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{smpdb-id}{Small Molecule Pathway Database identifier for this
 pathway.}
 \item{name}{Pathway name}
 \item{category}{Pathway category}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Metabolic, disease, and biological pathways that the drug is involved in, as
identified by the Small Molecule Protein Database (SMPDB).
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other pathway: 
\code{\link{drug_pathway_drugs}()},
\code{\link{drug_pathway_enzyme}()}
}
\concept{pathway}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_exp_prop}
\alias{drug_exp_prop}
\title{Drug Experimental Properties parser}
\usage{
drug_exp_prop(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{kind}{Name of the property.}
 \item{value}{Drug properties that have been experimentally proven.}
 \item{source}{Reference to the source of this experimental data.}
 \item{\emph{drugbank_id}}{drugbank id}
}

The following experimental properties are provided:
\describe{
 \item{Water Solubility}{The experimentally determined aqueous solubility
 of the molecule.}
 \item{Molecular Formula}{Protein formula of Biotech drugs}
 \item{Molecular Weight}{Protein weight of Biotech drugs.}
 \item{Melting Point}{The experimentally determined temperature at which the
  drug molecule changes from solid to liquid at atmospheric temperature.}
 \item{Boiling Point}{The experimentally determined temperature at which the
  drug molecule changes from liquid to gas at atmospheric temperature.}
 \item{Hydrophobicity}{The ability of a molecule to repel water rather than
 absorb or dissolve water.}
 \item{Isoelectric Point}{The pH value at which the net electric charge of a
 molecule is zero.}
 \item{caco2 Permeability}{A continuous line of heterogenous human epithelial
  colorectal adenocarcinoma cells, CAC02 cells are employed as a model of
  human intestinal absorption of various drugs and compounds. CAC02 cell
  permeability is ultimately an assay to measure drug absorption.}
 \item{pKa}{The experimentally determined pka value of the molecule}
 \item{logP}{The experimentally determined partition coefficient (LogP)
 based on the ratio of solubility of the molecule in 1-octanol compared to
 water.}
 \item{logS}{The intrinsic solubility of a given compound is the
 concentration in equilibrium with its solid phase that dissolves into
  solution, given as the natural logarithm (LogS) of the concentration.}
 \item{Radioactivity}{The property to spontaneously emit particles
 (alpha, beta, neutron) or radiation (gamma, K capture), or both at the same
 time, from the decay of certain nuclides.}
}
}
\description{
Drug properties that have been experimentally proven
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_groups}
\alias{drug_groups}
\title{Drug Groups parser}
\usage{
drug_groups(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 2 variables:
\describe{
 \item{group}{}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Groups that this drug belongs to. May include any of: approved, vet_approved,
 nutraceutical, illicit, withdrawn, investigational, and experimental.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbparser.R
\docType{package}
\name{dbparser}
\alias{dbparser}
\title{dbparser: A package for reading and parsing \strong{DrugBank} xml database.}
\description{
The main purpose of the `dbparser` package is to parse
[DrugBank](https://www.drugbank.ca/) database which is downloadable in XML format
from [this link](https://www.DrugBank.ca/releases/latest).
}
\details{
The parsed data can then be explored and analyzed as desired by the user
 with the ability to save parsed data into desired database as well.


To achieve this purpose, `dbparser`` package provides three main categories
 of functions:

- xml db reader,

- \strong{DrugBank} elements parsers,

- and database related methods.

For more information kindly check the
reference/index (https://docs.ropensci.org/dbparser/reference/index.html)
}
\section{xml db reader functions}{

 Reads \strong{DrugBank} xml database and build drug elements full tree in
 memory
}

\section{parsers functions}{

 Each parser function is responsible of parsing certain drug element and
 returning its tibble with the ability to save it in a predefined database.

 Check this tutorial
 (https://docs.ropensci.org/dbparser/articles/dbparser.html)
}

\section{database functions}{

 To open a connection to given database in order to store parsed
 \strong{DrugBank} elements database.

 Check this tutorial
 (https://docs.ropensci.org/dbparser/articles/Database_Saving.html)
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_manufacturer_node_parser.R
\name{drug_manufacturers}
\alias{drug_manufacturers}
\title{Drug Manufacturers parser}
\usage{
drug_manufacturers(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{generic}{A list of companies that are manufacturing the generic
  form of the drug.}
 \item{url}{A link to the companies that are manufacturing the drug.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
A list of companies that are manufacturing the commercially available forms
of this drug that are available in Canada and the Unites States.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_all_nodes_parser.R
\name{references}
\alias{references}
\title{Drugs/ Carriers/ Enzymes/ Targets/ Transporters references element parser}
\usage{
references(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\description{
Return a list of all references for drugs, carriers, enzymes, targets or
transporters
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other references: 
\code{\link{articles}},
\code{\link{attachments}},
\code{\link{books}},
\code{\link{links}}
}
\concept{references}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_intern_brand}
\alias{drug_intern_brand}
\title{Drug International Brands parser}
\usage{
drug_intern_brand(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 4 variables:
\describe{
 \item{brand}{The proprietary, well-known name for given to this drug by a
 manufacturer.}
 \item{company}{The company or manufacturer that uses this name.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
The proprietary names used by the manufacturers for commercially available
forms of the drug, focusing on brand names for products that are available
in countries other than Canada and the Unites States.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cett_polypeptide_others_parsers.R
\name{cett_poly_syn_doc}
\alias{cett_poly_syn_doc}
\alias{carriers_polypeptides_syn}
\alias{enzymes_polypeptides_syn}
\alias{targets_polypeptides_syn}
\alias{transporters_polypeptides_syn}
\title{Carriers/ Enzymes/ Targets/ Transporters Polypeptide Synonyms parsers}
\usage{
carriers_polypeptides_syn(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_polypeptides_syn(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_polypeptides_syn(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_polypeptides_syn(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 2 variables:
\describe{
  \item{synonym}{}
  \item{parent_key}{polypeptide id}
}
}
\description{
Extract descriptions of identified polypeptide synonyms for targets,
 enzymes, carriers, or transporters.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other cett: 
\code{\link{cett_actions_doc}},
\code{\link{cett_doc}},
\code{\link{cett_ex_identity_doc}},
\code{\link{cett_go_doc}},
\code{\link{cett_poly_doc}},
\code{\link{cett_poly_pfms_doc}}
}
\concept{cett}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_salts}
\alias{drug_salts}
\title{Drug Salts parser}
\usage{
drug_salts(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 1 variables:
\describe{
 \item{drugbank-id}{DrugBank identfiers of the available salt form(s).}
 \item{name}{Name of the available salt form(s)}
 \item{unii}{Unique Ingredient Identifier (UNII) of the available salt
 form(s).}
 \item{cas-number}{Chemical Abstracts Service (CAS) registry number assigned
  to the salt form(s) of the drug.}
 \item{inchikey}{IUPAC International Chemical Identifier (InChi) key
 identfier for the available salt form(s).}
 \item{average-mass}{Average molecular mass: the weighted average of the
  isotopic masses of the salt.}
 \item{monoisotopic-mass}{The mass of the most abundant isotope of the salt}
 \item{smiles}{The simplified molecular-input line-entry system (SMILES) is
 a line notation used for describing the structure of chemical species using
  short ASCII strings; calculated by ChemAxon.}
 \item{inchi}{A prediction of the IUPAC
 International Chemical Identifier (InChI); imported by ChemAxon.}
 \item{formula}{Indicates the simple numbers of each type of atom within the
  molecule; calculated by ChemAxon.}
 \item{\emph{drugbank_id}}{parent drugbank id}
}
}
\description{
Available salt forms of the drug. Ions such as hydrochloride, sodium,
 and sulfate are often added to the drug molecule to increase solubility,
 dissolution, or absorption.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/references_parsers.R
\name{attachments}
\alias{attachments}
\alias{drugs_attachments}
\alias{carriers_attachments}
\alias{enzymes_attachments}
\alias{targets_attachments}
\alias{transporters_attachments}
\title{Drugs/ Carriers/ Enzymes/ Targets/ Transporters attachments element parser}
\usage{
drugs_attachments(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

carriers_attachments(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_attachments(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_attachments(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_attachments(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 4 variables:
\describe{
  \item{ref-id}{Identifier for the article being referenced.
  This is unique across all reference types (books, links, article,
  attachments).}
  \item{title}{The title of the attachment.}
  \item{url}{The url to download the attachment from.}
  \item{\emph{parent_id}}{drug/carrier/target/enzyme/transporter id}
}
}
\description{
Return a list of attachment that were used as references for drugs carriers
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other references: 
\code{\link{articles}},
\code{\link{books}},
\code{\link{links}},
\code{\link{references}()}
}
\concept{references}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cett_polypeptide_others_parsers.R
\name{cett_go_doc}
\alias{cett_go_doc}
\alias{carriers_polypeptides_go}
\alias{enzymes_polypeptides_go}
\alias{targets_polypeptides_go}
\alias{transporters_polypeptides_go}
\title{Carriers/ Enzymes/ Targets/ Transporters Polypeptide GO Classifier
parsers}
\usage{
carriers_polypeptides_go(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_polypeptides_go(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_polypeptides_go(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_polypeptides_go(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 3 variables:
\describe{
  \item{category}{}
  \item{description}{}
  \item{parent_key}{polypeptide id}
}
}
\description{
Extract descriptions of identified polypeptide go classifier for targets,
 enzymes, carriers, or transporters.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other cett: 
\code{\link{cett_actions_doc}},
\code{\link{cett_doc}},
\code{\link{cett_ex_identity_doc}},
\code{\link{cett_poly_doc}},
\code{\link{cett_poly_pfms_doc}},
\code{\link{cett_poly_syn_doc}}
}
\concept{cett}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cett_polypeptide_general_infromation_parsers.R
\name{cett_poly_doc}
\alias{cett_poly_doc}
\alias{carriers_polypeptides}
\alias{enzymes_polypeptides}
\alias{targets_polypeptides}
\alias{transporters_polypeptides}
\title{Carriers/ Enzymes/ Targets/ Transporters Polypeptide parsers}
\usage{
carriers_polypeptides(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_polypeptides(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_polypeptides(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_polypeptides(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 20 variables:
\describe{
  \item{id}{\href{https://www.uniprot.org/}{Universal Protein Resource
  (UniProt) identifier}}
  \item{source}{Specifies whether the identified polypeptide ID is
  associated with any of the following UniProt knowledge bases:
  Swiss-Prot, which is manually annotated and reviewed, or TrEMBL,
  which is automatically annotated and not reviewed.}
  \item{name}{}
  \item{general_function}{General summary of the physiological function of
  the polypeptide}
  \item{specific_function}{A more specific description of the polypeptide’s
   physiological function within the cell.}
  \item{gene_name}{The short name commonly associated with the associated
   gene. Eg. PTGS1.}
  \item{locus}{The specific chromosomal location or position of the gene’s
   sequence on a chromosome.}
  \item{cellular_location}{The cellular location of the polypeptide.}
  \item{transmembrane_regions}{Areas of the polypeptide sequence that span
   a biological membrane.}
  \item{signal_regions}{Location of any signal peptides within the
   polypeptide sequence.}
  \item{theoretical_pi}{Theoretical isoelectric point.}
  \item{molecular_weight}{The molecular weight of the polypeptide.}
  \item{chromosome_location}{The chromosomal location of the polypeptide
  gene}
  \item{organism}{The organism in which this polypeptide functions.}
  \item{organism_ncbi_taxonomy_id}{}
  \item{amino_acid_sequence}{The amino acid sequence of the polypeptide}
  \item{amino_acid_format}{}
  \item{gene_sequence}{The sequence of the associated gene.}
  \item{gene_format}{}
  \item{parent_key}{carrier/ target/ enzyme/ transporter id}
}
}
\description{
Extract descriptions of identified polypeptide targets, enzymes, carriers,
 or transporters.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other cett: 
\code{\link{cett_actions_doc}},
\code{\link{cett_doc}},
\code{\link{cett_ex_identity_doc}},
\code{\link{cett_go_doc}},
\code{\link{cett_poly_pfms_doc}},
\code{\link{cett_poly_syn_doc}}
}
\concept{cett}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_snp_adverse_reactions}
\alias{drug_snp_adverse_reactions}
\title{Drug SNP Adverse Drug Reactions parser}
\usage{
drug_snp_adverse_reactions(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{protein-name}{Proteins involved in this SNP.}
 \item{gene-symbol}{Genes involved in this SNP.}
 \item{uniprot-id}{Universal Protein Resource (UniProt) identifiers for
 proteins involved in this pathway.}
 \item{rs-id}{The SNP Database identifier for this single nucleotide
  polymorphism.}
 \item{allele}{The alleles associated with the identified SNP.}
 \item{adverse-reaction}{}
 \item{description}{}
 \item{pubmed-id	}{Reference to PubMed article.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
The adverse drug reactions that may occur as a result of the listed single
nucleotide polymorphisms (SNPs)
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_price_node_parser.R
\name{drug_prices}
\alias{drug_prices}
\title{Drug Prices Parsers}
\usage{
drug_prices(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 5 variables:
\describe{
  \item{description}{}
  \item{cost}{Drug price per unit}
  \item{currency}{Currency of price, example: US.}
  \item{unit}{}
  \item{parent_id}{drugbank id}
}
}
\description{
Unit drug prices
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_all_nodes_parser.R
\name{cett}
\alias{cett}
\title{Run all CETT related parsers}
\usage{
cett(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a list of all drugs parsed tibbles
}
\description{
Run all parsers that retrieve carriers, enzymes, targets and transporters
 related information
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other collective_parsers: 
\code{\link{drugs}()},
\code{\link{run_all_parsers}()}
}
\concept{collective_parsers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_reaction_node_parser.R
\name{drug_reactions}
\alias{drug_reactions}
\title{Drug Reactions Parsers}
\usage{
drug_reactions(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 5 variables:
\describe{
  \item{sequence}{	Reactions are displayed within a numerical sequence}
  \item{left_drugbank_name}{The substrate of the reaction. Maybe a drug or a
   metabolite.}
  \item{rightt_drugbank_name}{	The product of the reaction. Maybe a drug or a
   metabolite.}
  \item{left_drugbank_id}{}
  \item{right_drugbank_id}{}
  \item{parent_id}{drugbank id}
}
}
\description{
Extract the sequential representation of the metabolic reactions that this
 drug molecule is involved in. Depending on available information, this may
 include metabolizing enzymes, reaction type, substrates, products,
 pharmacological activity of metabolites, and a structural representation of
  the biochemical reactions.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/references_parsers.R
\name{links}
\alias{links}
\alias{drugs_links}
\alias{carriers_links}
\alias{enzymes_links}
\alias{targets_links}
\alias{transporters_links}
\title{Drugs/ Carriers/ Enzymes/ Targets/ Transporters links element parser}
\usage{
drugs_links(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

carriers_links(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_links(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_links(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_links(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 4 variables:
\describe{
  \item{ref-id}{Name of the source website}
  \item{title}{Identifier for this drug in the given resource}
  \item{url}{The url of the website}
  \item{\emph{parent_id}}{drug/ carrier/ target/ enzyme/ transporter id}
}
}
\description{
Return a list of websites that were used as references for
Drugs/ Carriers/ Enzymes/ Targets/ Transporters
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other references: 
\code{\link{articles}},
\code{\link{attachments}},
\code{\link{books}},
\code{\link{references}()}
}
\concept{references}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_synonym_node_parser.R
\name{drug_syn}
\alias{drug_syn}
\title{Drug Synonyms parser}
\usage{
drug_syn(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 3 variables:
\describe{
 \item{language}{Names of the drug in languages other than English.}
 \item{coder}{Organisation or source providing the synonym. For example,
  INN indicates the synonym is an International Nonproprietary Name,
  while IUPAC indicates the synonym is the nomenclature designated by the
  International Union of Pure and Applied Chemistry.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Other names or identifiers that are associated with this drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_products}
\alias{drug_products}
\title{Drug Products parser}
\usage{
drug_products(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 32 variables:
\describe{
 \item{name}{The proprietary name(s) provided by the manufacturer for any
 commercially available products containing this drug.}
 \item{labeller}{The corporation responsible for labelling this product.}
 \item{ndc-id}{The National Drug Code (NDC) identifier of the drug}
 \item{ndc-product-code}{The National Drug Code (NDC) product code from the
  FDA National Drug Code directory.}
 \item{dpd-id}{Drug Product Database (DPD) identification number (a.k.a. DIN)
  from the Canadian Drug Product Database. Only present for drugs that are
  marketed in Canada}
 \item{ema-product-code}{EMA product code from the European Medicines Agency
 Database. Only present for products that are authorised by central procedure
  for marketing in the European Union.}
 \item{ema-ma-number}{EMA marketing authorisation number from the European
 Medicines Agency Database. Only present for products that are authorised by
  central procedure for marketing in the European Union.}
 \item{started-marketing-on}{The starting date for market approval.}
 \item{ended-marketing-on}{The ending date for market approval.}
 \item{dosage-form	}{The pharmaceutical formulation by which the drug is
 introduced into the body.}
 \item{strength}{The amount of active drug ingredient provided in the dosage}
 \item{route}{The path by which the drug or product is taken into the body}
 \item{fda-application-number}{The New Drug Application [NDA] number
 assigned to this drug by the FDA.}
 \item{over-the-counter}{A list of Over The Counter (OTC) forms of the drug.}
 \item{generic}{Whether this product is a generic drug.}
 \item{approved}{Indicates whether this drug has been approved by the
 regulating government.}
 \item{country}{The country where this commercially available drug has been
 approved.}
 \item{source}{Source of this product information. For example, a value of
 DPD indicates this information was retrieved from the Canadian Drug Product
  Database.}
 \item{standing}{One of good, discordant, or deprecated. Distinguishes
 products with up to date ingredient information (good) from products with
 conflicting information (discordant) or products that have been removed from
  an active label (deprecated).}
 \item{standing-updated-on}{The date on which the standing was last updated}
 \item{standing-reason}{Explains the non-good standing of the product.
 One of: ingredient_change, code_duplication, invalid, or removed.}
 \item{jurisdiction-marketing-category	}{The marketing category of this
 product in its jurisdiction}
 \item{branded}{Whether this product has a named brand}
 \item{prescription}{Whether this product is only available with
 a prescription}
 \item{unapproved}{Whether this product is not approved in its jurisdiction}
 \item{vaccine}{Whether this product is a vaccine}
 \item{allergenic}{Whether this product is used in allergenic testing}
 \item{cosmetic}{Whether this product is a cosmetic, such as sunscreen}
 \item{kit}{Whether this product is a kit composed of multiple distinct
 parts}
 \item{solo}{Whether this product has only a single active ingredient}
 \item{available}{Whether this product can be sold in its jurisdiction}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
A list of commercially available products in Canada and the United States
 that contain the drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_interactions}
\alias{drug_interactions}
\title{Drug Interactions parser}
\usage{
drug_interactions(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{drugbank-id	}{Drugbank ID of the interacting drug.}
 \item{name}{Name of the interacting drug.}
 \item{description}{Textual description of the physiological consequences
 of the drug interaction}
 \item{\emph{drugbank_id}}{parent drugbank id}
}
}
\description{
Drug-drug interactions detailing drugs that, when administered concomitantly
with the drug of interest, will affect its activity or result in adverse
effects. These interactions may be synergistic or antagonistic depending on
the physiological effects and mechanism of action of each drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_patents}
\alias{drug_patents}
\title{Drug Patents parser
A property right issued by the United States Patent and Trademark
Office (USPTO) to an inventor for a limited time, in exchange for public
disclosure of the invention when the patent is granted. Drugs may be issued
multiple patents.}
\usage{
drug_patents(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{number}{The patent number(s) associated with the drug.}
 \item{country}{The country that issued the patent rights.}
 \item{approved}{The date that the patent request was filed.}
 \item{expires}{The date that the patent rights expire.}
 \item{pediatric-extension}{Indicates whether or not a pediatric extension has been approved for
  the patent. Granted pediatric extensions provide an additional 6 months of
  market protection.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Drug Patents parser
A property right issued by the United States Patent and Trademark
Office (USPTO) to an inventor for a limited time, in exchange for public
disclosure of the invention when the patent is granted. Drugs may be issued
multiple patents.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_calc_prop}
\alias{drug_calc_prop}
\title{Drug Calculated Properties parser}
\usage{
drug_calc_prop(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 4 variables:
\describe{
 \item{kind}{Name of the property.}
 \item{value}{Predicted physicochemical properties; obtained by the use of
 prediction software such as ALGOPS and ChemAxon.}
 \item{source}{Name of the software used to calculate this property,
 either ChemAxon or ALOGPS.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Drug properties that have been predicted by ChemAxon or ALOGPS based on the
inputed chemical structure. Associated links below will redirect to
descriptions of the specific term.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_pathway_node_parsers.R
\name{drug_pathway_enzyme}
\alias{drug_pathway_enzyme}
\title{Drug Pathway Enzymes parser}
\usage{
drug_pathway_enzyme(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with pathway properties
}
\description{
Enzymes involved in this pathway.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other pathway: 
\code{\link{drug_pathway_drugs}()},
\code{\link{drug_pathway}()}
}
\concept{pathway}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_all_nodes_parser.R
\name{drugs}
\alias{drugs}
\title{Run all drug  related parsers}
\usage{
drugs(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a list of all drugs parsed tibbles
}
\description{
Run all parsers that retrieve drugs related information
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other collective_parsers: 
\code{\link{cett}()},
\code{\link{run_all_parsers}()}
}
\concept{collective_parsers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cett_general_information_parsers.R
\name{cett_doc}
\alias{cett_doc}
\alias{carriers}
\alias{enzymes}
\alias{targets}
\alias{transporters}
\title{Carriers/ Enzymes/ Targets/ Transporters parsers}
\usage{
carriers(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 6 variables (8 for enzymes):
\describe{
  \item{id}{Universal Protein Resource (UniProt) Identifier for the record}
  \item{name}{related name}
  \item{organism}{Organism that the protein comes from.}
  \item{known_action}{Whether the pharmacological action of the drug is due
   to this target interaction.}
  \item{inhibition-strength}{Whether the strength of enzyme inhibition is
  strong, moderate, or unknown. \strong{Only applies to enzymes}}
  \item{induction-strength}{Whether the strength of enzyme induction is
  strong or unknown. \strong{Only applies to enzymes}}
  \item{position}{related position}
  \item{parent_id}{drugbank id}
}
}
\description{
Protein targets of drug action, enzymes that are inhibited/induced or
involved in metabolism, and carrier or transporter proteins involved in
movement of the drug across biological membranes.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other cett: 
\code{\link{cett_actions_doc}},
\code{\link{cett_ex_identity_doc}},
\code{\link{cett_go_doc}},
\code{\link{cett_poly_doc}},
\code{\link{cett_poly_pfms_doc}},
\code{\link{cett_poly_syn_doc}}
}
\concept{cett}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_all_nodes_parser.R
\name{run_all_parsers}
\alias{run_all_parsers}
\title{extracts the all drug elements and return data as list of tibbles.}
\usage{
run_all_parsers(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
all drug elements tibbles
}
\description{
this functions extracts all element of drug nodes in \strong{DrugBank}
xml database with the option to save it in a predefined database via
passed database connection. it takes two optional arguments to
save the returned tibble in the database \code{save_table} and
\code{database_connection}.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other common: 
\code{\link{drug_element_options}()},
\code{\link{drug_element}()}

Other collective_parsers: 
\code{\link{cett}()},
\code{\link{drugs}()}
}
\concept{collective_parsers}
\concept{common}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_classfication_node_parser.R
\name{drug_classification}
\alias{drug_classification}
\title{Drug Classification parser}
\usage{
drug_classification(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 9 variables:
\describe{
  \item{description}{}
  \item{direct-parent}{}
  \item{kingdom}{}
  \item{superclass}{}
  \item{class}{}
  \item{subclass}{}
  \item{alternative-parent}{One or more alternative parents}
  \item{substituent}{One or more substituents}
  \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
A description of the hierarchical chemical classification of the drug;
imported from ClassyFire.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_sequence_node_parser.R
\name{drug_sequences}
\alias{drug_sequences}
\title{Drug Sequences parser}
\usage{
drug_sequences(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{sequence}{a textual representation of the sequence}
 \item{format}{Currently, only the FASTA format is used}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
The amino acid sequence; provided if the drug is a peptide.
}
\details{
Describes peptide sequences of biotech drugs
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_common_utilities.R
\name{read_drugbank_xml_db}
\alias{read_drugbank_xml_db}
\title{Reads \strong{DrugBank} xml database and load it into memory.}
\usage{
read_drugbank_xml_db(drugbank_db_path)
}
\arguments{
\item{drugbank_db_path}{\strong{string}, full path for the
\strong{DrugBank} xml or zip file.}
}
\value{
\strong{TRUE} when the loading process into memory to be used by
parser methods is completed successfully and \strong{FALSE} otherwise.
}
\description{
\code{read_drugbank_xml_db} loads \strong{DrugBank} xml database full tree
into memory.
}
\details{
This functions reads \strong{DrugBank} xml database and load it into memory
 for later processing. Hence; this method \strong{must} be called before any
 other function in the package and it needs to be called one time only.

It takes one single mandatory argument which is the location of DrugBank db.
}
\examples{
\dontrun{
read_drugbank_xml_db("db_full_path")
read_drugbank_xml_db(drugbank_db_path = "db_full_path")
}
}
\concept{DrugBank DB Loading}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cett_actions_parsers.R
\name{cett_actions_doc}
\alias{cett_actions_doc}
\alias{carriers_actions}
\alias{enzymes_actions}
\alias{targets_actions}
\alias{transporters_actions}
\title{Carriers/ Enzymes/ Targets/ Transporters Actions parsers}
\usage{
carriers_actions(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_actions(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_actions(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_actions(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 2 variables:
\describe{
  \item{action}{describe related action}
  \item{\emph{parent_id}}{carrier/ target/ enzyme/ transporter id}
}
}
\description{
Collection of related actions
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other cett: 
\code{\link{cett_doc}},
\code{\link{cett_ex_identity_doc}},
\code{\link{cett_go_doc}},
\code{\link{cett_poly_doc}},
\code{\link{cett_poly_pfms_doc}},
\code{\link{cett_poly_syn_doc}}
}
\concept{cett}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cett_polypeptide_others_parsers.R
\name{cett_poly_pfms_doc}
\alias{cett_poly_pfms_doc}
\alias{carriers_polypeptides_pfams}
\alias{enzymes_polypeptides_pfams}
\alias{targets_polypeptides_pfams}
\alias{transporters_polypeptides_pfams}
\title{Carriers/ Enzymes/ Targets/ Transporters Polypeptide PFAMS parsers}
\usage{
carriers_polypeptides_pfams(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_polypeptides_pfams(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_polypeptides_pfams(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_polypeptides_pfams(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 3 variables:
\describe{
  \item{name}{The sequence of the associated gene.}
  \item{identifier}{}
  \item{parent_key}{polypeptide id}
}
}
\description{
Extract descriptions of identified polypeptide PFAMS targets, enzymes,
 carriers, or transporters.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other cett: 
\code{\link{cett_actions_doc}},
\code{\link{cett_doc}},
\code{\link{cett_ex_identity_doc}},
\code{\link{cett_go_doc}},
\code{\link{cett_poly_doc}},
\code{\link{cett_poly_syn_doc}}
}
\concept{cett}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drugbank_metadata.R
\name{get_drugbank_version}
\alias{get_drugbank_version}
\title{Return uploaded drugbank database version}
\usage{
get_drugbank_version()
}
\value{
drugbank version
}
\description{
\code{get_drugbank_version} returns uploaded drugbank database version.
}
\examples{
\dontrun{
get_drugbank_version()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_all_nodes_parser.R
\name{drug_element_options}
\alias{drug_element_options}
\title{returns \code{drug_element} valid options.}
\usage{
drug_element_options()
}
\value{
list of \code{drug_element} valid options
}
\description{
returns \code{drug_element} valid options.
}
\examples{
\dontrun{
drug_element_options()
}
}
\seealso{
Other common: 
\code{\link{drug_element}()},
\code{\link{run_all_parsers}()}
}
\concept{common}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drugbank_metadata.R
\name{get_drugbank_exported_date}
\alias{get_drugbank_exported_date}
\title{Return uploaded drugbank database exported date}
\usage{
get_drugbank_exported_date()
}
\value{
drugbank exported date
}
\description{
\code{get_drugbank_exported_date} returns uploaded drugbank database
exported date.
}
\examples{
\dontrun{
get_drugbank_exported_date()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_external_links}
\alias{drug_external_links}
\title{Drug External Links parser}
\usage{
drug_external_links(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{resource}{Name of the source website.}
 \item{identifier}{Identifier for this drug in the given resource}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Links to other websites or databases providing information about this drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_ex_identity}
\alias{drug_ex_identity}
\title{Drug External Identifiers parser}
\usage{
drug_ex_identity(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{resource}{Name of the source database.}
 \item{identifier}{Identifier for this drug in the given resource.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Identifiers used in other websites or databases providing information about
this drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_atc_codes_node_parser.R
\name{drug_atc_codes}
\alias{drug_atc_codes}
\title{Drug ATC Codes element parser}
\usage{
drug_atc_codes(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 10 variables
}
\description{
The Anatomical Therapeutic Classification (ATC) code for the drug assigned
by the World Health Organization Anatomical Chemical Classification System.
}
\details{
Each `atc-code`` row has one or more level. The atc-code and level>
have a code  the code assigned by the World Health Organization Anatomical
Therapeutic Chemical Classification system.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_dosages}
\alias{drug_dosages}
\title{Drug Dosages parser}
\usage{
drug_dosages(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{form}{The pharmaceutical formulation by which the drug is introduced
 into the body}
 \item{route}{The path by which the drug or product is taken into the body.}
 \item{strength}{The amount of active drug ingredient provided in the dosage}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
A list of the commercially available dosages of the drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_pharmacology_parser.R
\name{drug_pharmacology}
\alias{drug_pharmacology}
\title{Drug Pharmacology parser}
\usage{
drug_pharmacology(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{indication}{The approved conditions, diseases, or states for which a
 drug can safely and effectively be used. An indication is considered to be
 FDA-approved when it has any of the following designations: NDA, ANDA, BLA,
  or OTC. May also include indications in other countries, such as Canada
  (through Health Canada) or in Europe
  (through the European Medicines Agency).}
 \item{pharmacodynamics}{A description of how the drug modifies or affects
 the organism it is being used in. May include effects in the body that are
  desired (enzyme or protein targets for example) and undesired
  (also known as “side effects”). This is in contrast to pharmacokinetics,
   which describes how the body modifies the drug being used.}
 \item{mechanism_of_action}{A component of pharmacodynamics that describes
 the biochemical interaction through which a drug produces its intended
 effect. May include the exact molecular protein or enzyme targets and/or
 a description of the physiological effects produced.}
 \item{toxicity}{Any adverse reaction, or side effect, that may or may not
 occur with use of the drug. May be attributed to a number of effects
 including: an enhanced therapeutic effect, rare anaphylactic reactions,
  interactions with other medications, or unanticipated binding of the
  molecule at different sites within the body.}
 \item{metabolism}{A description of the chemical degradation of the drug
 molecule within the body; most commonly by enzymes from the
 Cytochrome P450 (CYP) system in the liver.}
 \item{absorption}{A description of the movement of the drug from the site
  of administration into the bloodstream or target tissue. Common
  pharmacokinetic metrics used to evaluate absorption include Area Under
  the Curve (AUC), bioavailability (F), maximum concentration (Cmax), and
  time to maximum concentration (Tmax).}
 \item{half-life}{The period of time it takes for the amount of drug in the
 body to be reduced by one half. Provides a description of how quickly the
 drug is being eliminated and how much is available in the bloodstream.}
 \item{protein-binding	}{A description of the drug’s affinity for plama
 proteins and the proportion of the drug that is bound to them when in
 circulation within the body.}
 \item{route_of_elimination}{A description of the pathway that is used to
 excrete the drug from the body. Common pharmacokinetic parameters used to
 evaluate excretion include elemination half life, renal clearance, and
 tracking of radiolabelled compounds through the renal and GI system.}
 \item{volume_of_distribution}{The Vd of a drug represents the degree to
 which it is distributed into body tissue compared to the plasma.}
 \item{clearance}{A pharmacokinetic measurement of the rate of removal of the
  drug from plasma, expressed as mL/min; reflects the rate of elimination of
   the drug.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Describes the use, mechanism of action, pharmacokinetics, pharmacodynamics,
 and physiological or biochemical effects in the body.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/references_parsers.R
\name{books}
\alias{books}
\alias{drugs_textbooks}
\alias{carriers_textbooks}
\alias{enzymes_textbooks}
\alias{targets_textbooks}
\alias{transporters_textbooks}
\title{Drugs/ Carriers/ Enzymes/ Targets/ Transporters books element parser}
\usage{
drugs_textbooks(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

carriers_textbooks(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_textbooks(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_textbooks(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_textbooks(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 4 variables:
\describe{
  \item{ref-id}{Identifier for the article being referenced.
  This is unique across all reference types (books, links, article,
  attachments).}
  \item{isbn}{ISBN identifying the textbook.}
  \item{citation}{A Textbook citation in a standard format.}
  \item{\emph{parent_id}}{drug/ carrier/ target/ enzyme/ transporter id}
}
}
\description{
Return a list of text books that were used as references for drugs, carriers,
 enzymes, targets or transporters
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other references: 
\code{\link{articles}},
\code{\link{attachments}},
\code{\link{links}},
\code{\link{references}()}
}
\concept{references}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drugbank_metadata.R
\name{get_drugbank_metadata}
\alias{get_drugbank_metadata}
\title{Return uploaded drugbank database metadata}
\usage{
get_drugbank_metadata()
}
\value{
drugbank metadata
}
\description{
\code{get_drugbank_metadata} returns uploaded drugbank database version and
exported date.
}
\examples{
\dontrun{
get_drugbank_metadata()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_pdb_entries}
\alias{drug_pdb_entries}
\title{Drug pdb-entries parser}
\usage{
drug_pdb_entries(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{pdb-entry}{}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Protein Data Bank (PDB) identifiers for this drug.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_affected_organisms}
\alias{drug_affected_organisms}
\title{Drug Affected Organism parser}
\usage{
drug_affected_organisms(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 2 variables:
\describe{
 \item{affected-organism}{affected-organism name}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
Organisms in which the drug may display activity; activity may depend on
local susceptibility patterns and resistance.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_main_node_parser.R
\name{drug_general_information}
\alias{drug_general_information}
\title{Drugs General Information parser}
\usage{
drug_general_information(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 15 variables:
\describe{
  \item{primary_key}{Drugbank id}
  \item{other_keys}{Other identifiers that may be associated with the drug}
  \item{type}{	Either small molecule, or biotech. Biotech is used for any
  drug that is derived from living systems or organisms, usually composed of
   high molecular weight mixtures of protein, while small molecule describes
    a low molecular weight organic compound.}
  \item{name}{}
  \item{created}{Date that this drug was first added to DrugBank.}
  \item{updated}{Denotes when this drug was last updated in DrugBank.}
  \item{description}{Descriptions of drug chemical properties,
   history and regulatory status.}
  \item{cas_number}{The Chemical Abstracts Service (CAS) registry number
   assigned to the drug.}
  \item{\emph{unii}}{Unique Ingredient Identifier (UNII) of this drug.}
  \item{average_mass}{The weighted average of the isotopic masses of the
  drug}
  \item{state}{One of solid, liquid, or gas}
  \item{monoisotopic_mass}{The mass of the most abundant isotope of the drug}
  \item{synthesis_reference}{Citation for synthesis of the drug molecule.}
  \item{fda_label}{Contains a URL for accessing the uploaded United States Food
  and Drug Administration (FDA) Monograph for this drug.}
  \item{msds}{Contains a URL for accessing the Material Safety Data Sheet
  (MSDS) for this drug.}
}
}
\description{
A description of the hierarchical chemical classification of the drug;
imported from ClassyFire.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cett_polypeptide_others_parsers.R
\name{cett_ex_identity_doc}
\alias{cett_ex_identity_doc}
\alias{carriers_polypep_ex_ident}
\alias{enzymes_polypep_ex_ident}
\alias{targets_polypep_ex_ident}
\alias{transporters_polypep_ex_ident}
\title{Carriers/ Enzymes/ Targets/ Transporters Polypeptide External Identifiers
parsers}
\usage{
carriers_polypep_ex_ident(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_polypep_ex_ident(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_polypep_ex_ident(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_polypep_ex_ident(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 3 variables:
\describe{
  \item{resource}{Name of the source database.}
  \item{identifier}{Identifier for this drug in the given resource.}
  \item{parent_key}{polypeptide id}
}
}
\description{
Extract descriptions of identified polypeptide external identifiers for
targets, enzymes, carriers, or transporters.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other cett: 
\code{\link{cett_actions_doc}},
\code{\link{cett_doc}},
\code{\link{cett_go_doc}},
\code{\link{cett_poly_doc}},
\code{\link{cett_poly_pfms_doc}},
\code{\link{cett_poly_syn_doc}}
}
\concept{cett}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/references_parsers.R
\name{articles}
\alias{articles}
\alias{drugs_articles}
\alias{carriers_articles}
\alias{enzymes_articles}
\alias{targets_articles}
\alias{transporters_articles}
\title{Drugs/ Carriers/ Enzymes/ Targets/ Transporters articles element parser}
\usage{
drugs_articles(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

carriers_articles(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

enzymes_articles(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

targets_articles(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)

transporters_articles(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 4 variables:
\describe{
  \item{ref-id}{Identifier for the article being referenced.
  This is unique across all reference types (books, links, article,
  attachments).}
  \item{pubmed-id}{The PubMed identifier for the article.}
  \item{citation}{Article citation in a standard format.}
  \item{\emph{parent_id}}{drug/carrier/target/enzyme/transporter id}
}
}
\description{
Return a list of articles that were used as references for drugs carriers
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other references: 
\code{\link{attachments}},
\code{\link{books}},
\code{\link{links}},
\code{\link{references}()}
}
\concept{references}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_all_nodes_parser.R
\name{drug_element}
\alias{drug_element}
\title{extracts the given drug elements and return data as list of tibbles.}
\usage{
drug_element(
  elements_options = c("all"),
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{elements_options}{list,  options of elements to be parsed. default is
"all"}

\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.
@return list of selected drug elements tibbles}
}
\description{
\code{drug_element} returns list of tibbles of drugs selected
elements.
}
\details{
this functions extracts selected element of drug nodes in \strong{DrugBank}
xml database with the option to save it in a predefined database via
passed database connection. it takes two optional arguments to
save the returned tibble in the database \code{save_table} and
 \code{database_connection}.
it must be called after \code{\link{read_drugbank_xml_db}} function like
any other parser function.
if \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.

drug_element_options can be called to know the valid options for
this method
}
\examples{
\dontrun{
# return only the parsed tibble
drug_element()

# will throw an error, as database_connection is NULL
drug_element(save_table = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist before read it and return its data.
drug_element(save_csv = TRUE)


sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
drug_element(save_table = TRUE, database_connection = sqlite_con)

# save in database, save parsed tibble as csv if it does not
# exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
drug_element(save_table = TRUE, save_csv = TRUE,
 database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location and
# return parsed tibble.
# if the csv exist before read it and return its data.
drug_element(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current
# location and return parsed tibble.
# if the csv exist override it and return it.
drug_element(save_csv = TRUE, csv_path = TRUE, override = TRUE)
drug_element(c("drug_ahfs_codes", "drug_carriers"), save_table = TRUE)
drug_element(save_table = FALSE)
drug_element(c("drug_ahfs_codes", "drug_carriers"))
}
}
\seealso{
Other common: 
\code{\link{drug_element_options}()},
\code{\link{run_all_parsers}()}
}
\concept{common}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_packagers}
\alias{drug_packagers}
\title{Drug Packagers parser}
\usage{
drug_packagers(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 2 variables:
\describe{
 \item{name}{}
 \item{url}{A link to any companies that are packaging the drug for
 re-distribution.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
A list of companies that are packaging the drug for re-distribution.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_pathway_node_parsers.R
\name{drug_pathway_drugs}
\alias{drug_pathway_drugs}
\title{Drug Pathway Drugs parser}
\usage{
drug_pathway_drugs(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with pathway drugsproperties
}
\description{
Drugs involved in this pathway.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other pathway: 
\code{\link{drug_pathway_enzyme}()},
\code{\link{drug_pathway}()}
}
\concept{pathway}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_reaction_node_parser.R
\name{drug_reactions_enzymes}
\alias{drug_reactions_enzymes}
\title{Drug Reactions Enzymes Parsers}
\usage{
drug_reactions_enzymes(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with 3 variables:
\describe{
  \item{name}{}
  \item{uniprot-id}{}
  \item{parent_id}{drugbank id}
}
}
\description{
EEnzymes involved in metabolizing this drug
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_snp_effects}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drug_parsers.R
\name{drug_snp_effects}
\alias{drug_snp_effects}
\title{Drug SNP Effects parser}
\usage{
drug_snp_effects(
  save_table = FALSE,
  save_csv = FALSE,
  csv_path = ".",
  override_csv = FALSE,
  database_connection = NULL
)
}
\arguments{
\item{save_table}{boolean, save table in database if true.}

\item{save_csv}{boolean, save csv version of parsed tibble if true}

\item{csv_path}{location to save csv files into it, default is current
location, save_csv must be true}

\item{override_csv}{override existing csv, if any, in case it is true in the
new parse operation}

\item{database_connection}{DBI connection object that holds a connection to
user defined database. If \code{save_table} is enabled without providing
value for this function an error will be thrown.}
}
\value{
a tibble with the following variables:
\describe{
 \item{protein-name}{Proteins involved in this SNP.}
 \item{gene-symbol}{Genes involved in this SNP.}
 \item{uniprot-id}{Universal Protein Resource (UniProt) identifiers for
 proteins involved in this pathway.}
 \item{rs-id}{	The SNP Database identifier for this single nucleotide
 polymorphism.}
 \item{allele}{The alleles associated with the identified SNP.}
 \item{defining-change}{}
 \item{description}{A written description of the SNP effects.}
 \item{pubmed-id	}{Reference to PubMed article.}
 \item{\emph{drugbank_id}}{drugbank id}
}
}
\description{
A list of single nucleotide polymorphisms (SNPs) relevent to drug activity or
 metabolism, and the effects these may have on pharmacological activity.
 SNP effects in the patient may require close monitoring, an increase or
 decrease in dose, or a change in therapy.
}
\section{read_drugbank_xml_db}{

\code{\link{read_drugbank_xml_db}} function must be called first before any
parser.

If \code{\link{read_drugbank_xml_db}} is called before for any reason, so
no need to call it again before calling this function.
}

\examples{
\dontrun{
# the same parameters and usage will be applied for any parser
# return only the parsed tibble
run_all_parsers()

# will throw an error, as database_connection is NULL
run_all_parsers(save_table = TRUE)

# save in database in SQLite in memory database and return parsed tibble
sqlite_con <- DBI::dbConnect(RSQLite::SQLite())
run_all_parsers(save_table = TRUE, database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in current location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE)

# save in database, save parsed tibble as csv,
# if it does not exist in current location and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_table = TRUE, save_csv = TRUE,
database_connection = sqlite_con)

# save parsed tibble as csv if it does not exist in given location,
# and return parsed tibble.
# if the csv exist before read it and return its data.
run_all_parsers(save_csv = TRUE, csv_path = TRUE)

# save parsed tibble as csv if it does not exist in current location and
# return parsed tibble.
# if the csv exist override it and return it.
run_all_parsers(save_csv = TRUE, csv_path = TRUE, override = TRUE)
}
}
\seealso{
Other drugs: 
\code{\link{drug_affected_organisms}()},
\code{\link{drug_ahfs_codes}()},
\code{\link{drug_atc_codes}()},
\code{\link{drug_calc_prop}()},
\code{\link{drug_categories}()},
\code{\link{drug_classification}()},
\code{\link{drug_dosages}()},
\code{\link{drug_ex_identity}()},
\code{\link{drug_exp_prop}()},
\code{\link{drug_external_links}()},
\code{\link{drug_food_interactions}()},
\code{\link{drug_general_information}()},
\code{\link{drug_groups}()},
\code{\link{drug_interactions}()},
\code{\link{drug_intern_brand}()},
\code{\link{drug_manufacturers}()},
\code{\link{drug_mixtures}()},
\code{\link{drug_packagers}()},
\code{\link{drug_patents}()},
\code{\link{drug_pdb_entries}()},
\code{\link{drug_pharmacology}()},
\code{\link{drug_prices}()},
\code{\link{drug_products}()},
\code{\link{drug_reactions_enzymes}()},
\code{\link{drug_reactions}()},
\code{\link{drug_salts}()},
\code{\link{drug_sequences}()},
\code{\link{drug_snp_adverse_reactions}()},
\code{\link{drug_syn}()}
}
\concept{drugs}
