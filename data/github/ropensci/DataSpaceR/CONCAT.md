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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# DataSpaceR <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/DataSpaceR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/DataSpaceR/actions)
[![codecov](https://codecov.io/gh/ropensci/DataSpaceR/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/DataSpaceR/branch/main)
[![CRAN
Status](https://www.r-pkg.org/badges/version/DataSpaceR)](https://cran.r-project.org/package=DataSpaceR)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/DataSpaceR)](https://www.rpackages.io/package/DataSpaceR)
[![monthly](https://cranlogs.r-pkg.org/badges/DataSpaceR)](https://www.rpackages.io/package/DataSpaceR)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![](https://badges.ropensci.org/261_status.svg)](https://github.com/ropensci/software-review/issues/261)
<!-- badges: end -->

`DataSpaceR` is an R interface to [the CAVD
DataSpace](https://dataspace.cavd.org), a data sharing and discovery
tool that facilitates exploration of HIV immunological data from
pre-clinical and clinical HIV vaccine studies.

This package is intended for use by immunologists, bioinformaticians, and
statisticians in HIV vaccine research, or anyone interested in the
analysis of HIV immunological data across assays, studies, and time.

This package simplifies access to the database by taking advantage of
the standardization of the database to hide all the
[Rlabkey](https://cran.r-project.org/package=Rlabkey) specific code away
from the user, and it allows the users to access the study-specific
datasets via [an object-oriented
paradigm](https://cran.r-project.org/package=R6/readme/README.html).

## Examples & Documentation

For more detailed examples and detailed documentation, see [the
introductory
vignette](https://docs.ropensci.org/DataSpaceR/articles/DataSpaceR.html)
and [the pkgdown site](https://docs.ropensci.org/DataSpaceR/).

For a quick guide of how to use the API, see our [cheat sheet](https://dataspace.cavd.org/_webdav/static/%40files/documents/dataspacer_cheat_sheet.pdf) .

## Installation

Install from CRAN:

``` r
install.packages("DataSpaceR")
```

You can install the latest development version from
[GitHub](https://github.com/ropensci/DataSpaceR) with
[devtools](https://cran.r-project.org/package=devtools):

``` r
# install.packages("devtools")
devtools::install_github("ropensci/DataSpaceR")
```

## Register and set DataSpace credential

The database is accessed with the user’s credentials. A netrc file
storing login and password information is ***required***.

1.  [Create an account](https://dataspace.cavd.org/) and read the terms
    of use
2.  On your R console, create a netrc file using a function from
    `DataSpaceR`:

``` r
library(DataSpaceR)
writeNetrc(
  login = "yourEmail@address.com", 
  password = "yourSecretPassword",
  netrcFile = "/your/home/directory/.netrc" # use getNetrcPath() to get the default path 
)
```

This will create a netrc file in your home directory.

***Alternatively***, you can manually create a netrc file in the
computer running R.

-   On Windows, this file should be named `_netrc`
-   On UNIX, it should be named `.netrc`
-   The file should be located in the user’s home directory, and the
    permissions on the file should be unreadable for everybody except
    the owner
-   To determine home directory, run `Sys.getenv("HOME")` in R

The following three lines must be included in the `.netrc` or `_netrc`
file either separated by white space (spaces, tabs, or newlines) or
commas. Multiple such blocks can exist in one file.

    machine dataspace.cavd.org
    login myuser@domain.com
    password supersecretpassword

See
[here](https://www.labkey.org/wiki/home/Documentation/page.view?name=netrc)
for more information about `netrc`.

## Usage

The general idea is that the user:

1.  creates an instance of `DataSpaceConnection` class via `connectDS`
2.  browses available studies and groups in the instance via
    `availableStudies` and `availableGroups`
3.  creates a connection to a specific study via `getStudy` or a group
    via `getGroup`
4.  retrieves datasets by name via `getDataset`

### for example:

``` r
library(DataSpaceR)
#> By exporting data from the CAVD DataSpace, you agree to be bound by the Terms of Use available on the CAVD DataSpace sign-in page at https://dataspace.cavd.org

con <- connectDS()
con
#> <DataSpaceConnection>
#>   URL: https://dataspace.cavd.org
#>   User: jkim2345@scharp.org
#>   Available studies: 273
#>     - 77 studies with data
#>     - 5049 subjects
#>     - 423195 data points
#>   Available groups: 6
#>   Available publications: 1530
#>     - 12 publications with data
```

`connectDS()` will create a connection to DataSpace.

### available studies can be listed by `availableStudies` field

``` r
knitr::kable(head(con$availableStudies))
```

| study\_name | short\_name                    | title                                                                                                                                                                                                                                                              | type               | status   | stage            | species            | start\_date | strategy                             | network | data\_availability | ni\_data\_availability |
|:------------|:-------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------|:---------|:-----------------|:-------------------|:------------|:-------------------------------------|:--------|:-------------------|:-----------------------|
| cor01       | NA                             | The correlate of risk targeted intervention study (CORTIS): A randomized, partially-blinded, clinical trial of isoniazid and rifapentine (3HP) therapy to prevent pulmonary tuberculosis in high-risk individuals identified by a transcriptomic correlate of risk | Phase III          | Inactive | Assays Completed | Human              | NA          | NA                                   | GH-VAP  | NA                 | NA                     |
| cvd232      | Parks\_RV\_232                 | ​Limiting Dose Vaginal SIVmac239 Challenge of RhCMV-SIV vaccinated Indian rhesus macaques.                                                                                                                                                                         | Pre-Clinical NHP   | Inactive | Assays Completed | Rhesus macaque     | 2009-11-24  | Vector vaccines (viral or bacterial) | CAVD    | NA                 | NA                     |
| cvd234      | Zolla-Pazner\_Mab\_test1 Study | Zolla-Pazner\_Mab\_Test1                                                                                                                                                                                                                                           | Antibody Screening | Inactive | Assays Completed | Non-Organism Study | 2009-02-03  | Prophylactic neutralizing Ab         | CAVD    | NA                 | NA                     |
| cvd235      | mAbs potency                   | Weiss mAbs potency                                                                                                                                                                                                                                                 | Antibody Screening | Inactive | Assays Completed | Non-Organism Study | 2008-08-21  | Prophylactic neutralizing Ab         | CAVD    | NA                 | NA                     |
| cvd236      | neutralization assays          | neutralization assays                                                                                                                                                                                                                                              | Antibody Screening | Active   | In Progress      | Non-Organism Study | 2009-02-03  | Prophylactic neutralizing Ab         | CAVD    | NA                 | NA                     |
| cvd238      | Gallo\_PA\_238                 | HIV-1 neutralization responses in chronically infected individuals                                                                                                                                                                                                 | Antibody Screening | Inactive | Assays Completed | Non-Organism Study | 2009-01-08  | Prophylactic neutralizing Ab         | CAVD    | NA                 | NA                     |

### available groups can be listed by `availableGroups` field

``` r
knitr::kable(con$availableGroups)
```

| group\_id | label                              | original\_label                    | description                                                                                                               | created\_by | shared |   n | studies                        |
|----------:|:-----------------------------------|:-----------------------------------|:--------------------------------------------------------------------------------------------------------------------------|:------------|:-------|----:|:-------------------------------|
|       216 | mice                               | mice                               | NA                                                                                                                        | readjk      | FALSE  |  75 | cvd468, cvd483, cvd316, cvd331 |
|       217 | CAVD 242                           | CAVD 242                           | This is a fake group for CAVD 242                                                                                         | readjk      | FALSE  |  30 | cvd242                         |
|       220 | NYVAC durability comparison        | NYVAC\_durability                  | Compare durability in 4 NHP studies using NYVAC-C (vP2010) and NYVAC-KC-gp140 (ZM96) products.                            | ehenrich    | TRUE   |  78 | cvd281, cvd434, cvd259, cvd277 |
|       224 | cvd338                             | cvd338                             | NA                                                                                                                        | readjk      | FALSE  |  36 | cvd338                         |
|       228 | HVTN 505 case control subjects     | HVTN 505 case control subjects     | Participants from HVTN 505 included in the case-control analysis                                                          | drienna     | TRUE   | 189 | vtn505                         |
|       230 | HVTN 505 polyfunctionality vs BAMA | HVTN 505 polyfunctionality vs BAMA | Compares ICS polyfunctionality (CD8+, Any Env) to BAMA mfi-delta (single Env antigen) in the HVTN 505 case control cohort | drienna     | TRUE   | 170 | vtn505                         |

***Note***: A group is a curated collection of participants from
filtering of treatments, products, studies, or species, and it is
created in [the DataSpace
App](https://dataspace.cavd.org/cds/CAVD/app.view).

Check out [the reference
page](https://docs.ropensci.org/DataSpaceR/reference/DataSpaceConnection.html)
of `DataSpaceConnection` for all available fields and methods.

### create an instance of `cvd408`

``` r
cvd408 <- con$getStudy("cvd408")
cvd408
#> <DataSpaceStudy>
#>   Study: cvd408
#>   URL: https://dataspace.cavd.org/CAVD/cvd408
#>   Available datasets:
#>     - Binding Ab multiplex assay
#>     - Demographics
#>     - Intracellular Cytokine Staining
#>     - Neutralizing antibody
#>   Available non-integrated datasets:
class(cvd408)
#> [1] "DataSpaceStudy" "R6"
```

### available datasets can be listed by `availableDatasets` field

``` r
knitr::kable(cvd408$availableDatasets)
```

| name         | label                           |    n | integrated |
|:-------------|:--------------------------------|-----:|:-----------|
| BAMA         | Binding Ab multiplex assay      | 1080 | TRUE       |
| Demographics | Demographics                    |   20 | TRUE       |
| ICS          | Intracellular Cytokine Staining | 3720 | TRUE       |
| NAb          | Neutralizing antibody           |  540 | TRUE       |

which will print names of available datasets.

### Neutralizing Antibody dataset (`NAb`) can be retrieved by:

``` r
NAb <- cvd408$getDataset("NAb")
dim(NAb)
#> [1] 540  33
colnames(NAb)
#>  [1] "participant_id"      "participant_visit"   "visit_day"          
#>  [4] "assay_identifier"    "summary_level"       "specimen_type"      
#>  [7] "antigen"             "antigen_type"        "virus"              
#> [10] "virus_type"          "virus_insert_name"   "clade"              
#> [13] "neutralization_tier" "tier_clade_virus"    "target_cell"        
#> [16] "initial_dilution"    "titer_ic50"          "titer_ic80"         
#> [19] "response_call"       "nab_lab_source_key"  "lab_code"           
#> [22] "exp_assayid"         "titer_id50"          "titer_id80"         
#> [25] "nab_response_id50"   "nab_response_id80"   "slope"              
#> [28] "vaccine_matched"     "study_prot"          "virus_full_name"    
#> [31] "virus_species"       "virus_host_cell"     "virus_backbone"
```

Check out [the reference
page](https://docs.ropensci.org/DataSpaceR/reference/DataSpaceStudy.html)
of `DataSpaceStudy` for all available fields and methods.

***Note***: The package uses a
[R6](https://cran.r-project.org/package=R6) class to represent the
connection to a study and get around some of R’s copy-on-change
behavior.

## Meta

-   Please [report any issues or
    bugs](https://github.com/ropensci/DataSpaceR/issues).
-   License: GPL-3
-   Get citation information for `DataSpaceR` in R doing
    `citation(package = 'DataSpaceR')`
-   Please note that this project is released with a [Contributor Code
    of
    Conduct](https://github.com/ropensci/DataSpaceR/blob/main/CONDUCT.md).
    By participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# DataSpaceR 0.7.5

* Fixed bug concerning missing cookie session name for app reports which will update this package to work when run on CDS as in a report.
* Show NI data in available datasets active binding which will allow users to look for NI (non-integrated) data in DSR.
* Set active binding names to snake case which is in effort to standardize fields names used in the API.
* Add all authors to publication table which will allow users to find publication by authors who were not listed as the primary author.

# DataSpaceR 0.7.4

* Added `virusMetadata` field to `DataSpaceConnection` which shows virus metadata for viruses used in NAb assays (thanks @helenmiller16 #26)
* Added `availablePublications` field to `DataSpaceConnection` which summarizes all available publications  (thanks @helenmiller16 #27)
* Added `downloadPublicationData` method to `DataSpaceConnection` which will download available publication data for a specified publication (thanks @helenmiller16 #27)
* Added a vignette: `Accessing Publication Data` (thanks @helenmiller16 #27)
* Made `PKMAb` dataset available to retrieve in `DataSpaceStudy` (thanks @helenmiller16 and @jmtaylor-fhcrc #27)
* Relaxed some assumptions when pulling non-integrated data, which allows users to pull non-integrated mab data, like in cvd812. (thanks @helenmiller16 #28)
* Fixed a test for `mabGridSummary$geometric_mean_curve_ic50` calculation to reflect updated data. (thanks @helenmiller16 #28)
* Updated documentation using the R6 documentation syntax in roxygen

# DataSpaceR 0.7.3

* Remove `tolower` in functions that check study names. (thanks @jmtaylor-fhcrc #23)
* Update DataSpaceStudy methods for non-integrated datasets. (thanks @helenmiller16 #24)

# DataSpaceR 0.7.2

* Fix broken and invalid URLs for CRAN submission

# DataSpaceR 0.7.1

* Fix broken and invalid URLs for CRAN submission

# DataSpaceR 0.7.0

* Added `DataSpaceMab` class and several methods and fields in `DataSpaceConnection` to allow access to monoclonal antibody data (thanks @jmtaylor-fhcrc #14, #19, #21)

# DataSpaceR 0.6.3

* Prepare for CRAN submission

# DataSpaceR 0.6.2

* Write the netrc in temporary directory as default in `writeNetrc`. (CRAN requirement)

# DataSpaceR 0.6.1

* Remove LICENSE file for CRAN submission

# DataSpaceR 0.6.0

* rOpenSci submission: https://github.com/ropensci/software-review/issues/261

# DataSpaceR 0.5.2

* Modified `getDataset` method in `DataSpaceStudy` to take `mergeExtra` as an argument in order to merge extra information (demographics and treatment arm). #5

# DataSpaceR 0.5.1

* Adjusted the package to use the latest version of Rlabkey (v2.2) and httr packages.

# DataSpaceR 0.5.0

* Modified `DataSpaceConnection` and `DataSpaceStudy` classes to convert data.frame objects to data.table.
* Added a package startup message on the terms of use.
* Created `getGroup` method in `DataSpaceConnection` class and deprecated `groupId` in `getStudy` method.

# DataSpaceR 0.4.2

* Added `refresh` method for both connection and study classes.
* Added `studyInfo` field in study class.

# DataSpaceR 0.4.1

* Included additional columns (`short_name`, `type`, `status`, `stage`, `species`, `start_date`, `strategy`) in `availableStudies`.
* Updated the introductory vignette.

# DataSpaceR 0.4.0

* Added `availableGroups` field to `DataSpaceConnection` class.
* Modified `getStudy` method in `DataSpaceConnection` to take `groupId` as an argument in order to create a `DataSpaceStudy` object for a particular group.

# DataSpaceR 0.3.0

* Updated `DataSpaceConnection` class and `connectDS` constructor to be not study-specific.

# DataSpaceR 0.2.1

* Modified `getUserEmail()` to get email from netrc file.
* Updated license to GPL-3.

# DataSpaceR 0.2.0

* Added `getVariableInfo` method to `DataSpaceConnection` class.
* Changed variable names in `getAvailableDatasets` to lowercase.
* Added `treatmentArm` field to `DataSpaceConnection` class.
* Removed `getAvailableDatasets` method (now a private method).
* Changed the default connection from staging (`dataspace-staging.cavd.org`) to production (`dataspace.cavd.org`).
* Added an option to connect to the staging server.
* Renamed `write_netrc` and `check_netrc` to `writeNetrc` and `checkNetrc`.

# DataSpaceR 0.1.0

* Initialized the package with `connectDS()` and other basic functions.
* Added a test framework on `tests` to test the package.
* Added `.travis.yml` to build and check the package automatically. (not used yet)
* Added `codecov.yml` to track code coverage. (not used yet)
* Added `README.Rmd` to introduce the package.
* Added `NEWS.md` to track changes to the package.
* Added `_pkgdown.yml` and `/docs` to build a package website.
* Added a introductory vignette called `Intro_to_DataSpaceR.Rmd`.
## Release summary

* This is a new release.
* Jason Taylor (jmtaylor@fredhutch.org) is the new maintainer.

## Test environments
* local R installation, R 4.1.1
* ubuntu 18.04 and macOS (on GitHub Actions), R 4.1.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Jason Taylor <jmtaylor@fredhutch.org>'

  Jason Taylor <jmtaylor@fredhutch.org>
New maintainer:
Old maintainer(s):
  Ju Yeong Kim <jkim2345@fredhutch.org>
```
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

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

Please note that the DataSpaceR project is released with a
[Contributor Code of Conduct](CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://ropensci.github.io/dev_guide/contributingguide.html)
for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email? 

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md 
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
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

# DataSpaceR <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->
[![R build status](https://github.com/ropensci/DataSpaceR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/DataSpaceR/actions)
[![codecov](https://codecov.io/gh/ropensci/DataSpaceR/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/DataSpaceR/branch/main)
[![CRAN Status](https://www.r-pkg.org/badges/version/DataSpaceR)](https://cran.r-project.org/package=DataSpaceR)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/DataSpaceR)](https://www.rpackages.io/package/DataSpaceR)
[![monthly](https://cranlogs.r-pkg.org/badges/DataSpaceR)](https://www.rpackages.io/package/DataSpaceR)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![](https://badges.ropensci.org/261_status.svg)](https://github.com/ropensci/software-review/issues/261)
<!-- badges: end -->



`DataSpaceR` is an R interface to [the CAVD DataSpace](https://dataspace.cavd.org), a data sharing and discovery tool that facilitates exploration of HIV immunological data from pre-clinical and clinical HIV vaccine studies.

The package is intended for use by immunologists, bioinformaticians, and statisticians in HIV vaccine research, or anyone interested in the analysis of HIV immunological data across assays, studies, and time.

This package simplifies access to the database by taking advantage of the standardization of the database to hide all the [Rlabkey](https://cran.r-project.org/package=Rlabkey) specific code away from the user, and it allows the users to access the study-specific datasets via [an object-oriented paradigm](https://cran.r-project.org/package=R6/readme/README.html).


## Examples & Documentation

For more detailed examples and detailed documentation, see [the introductory vignette](https://docs.ropensci.org/DataSpaceR/articles/DataSpaceR.html) and [the pkgdown site](https://docs.ropensci.org/DataSpaceR/).


## Installation

Install from CRAN:

```{r, eval=FALSE}
install.packages("DataSpaceR")
```

You can install the latest development version from [GitHub](https://github.com/ropensci/DataSpaceR) with [devtools](https://cran.r-project.org/package=devtools):

```{r, eval=FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/DataSpaceR")
```


## Register and set DataSpace credential

The database is accessed with the user's credentials. A netrc file storing 
login and password information is ***required***.

1. [Create an account](https://dataspace.cavd.org/) and read the terms of use
1. On your R console, create a netrc file using a function from `DataSpaceR`:

```{r, eval=FALSE}
library(DataSpaceR)
writeNetrc(
  login = "yourEmail@address.com", 
  password = "yourSecretPassword",
  netrcFile = "/your/home/directory/.netrc" # use getNetrcPath() to get the default path 
)
```

This will create a netrc file in your home directory.

***Alternatively***, you can manually create a netrc file in the computer running R.

* On Windows, this file should be named `_netrc`
* On UNIX, it should be named `.netrc`
* The file should be located in the user's home directory, and the permissions on the file should be unreadable for everybody except the owner
* To determine home directory, run `Sys.getenv("HOME")` in R

The following three lines must be included in the `.netrc` or `_netrc` file either separated by white space (spaces, tabs, or newlines) or commas. Multiple such blocks can exist in one file.

```
machine dataspace.cavd.org
login myuser@domain.com
password supersecretpassword
```

See [here](https://www.labkey.org/wiki/home/Documentation/page.view?name=netrc) for more information about `netrc`.


## Usage

The general idea is that the user:

1. creates an instance of `DataSpaceConnection` class via `connectDS`
1. browses available studies and groups in the instance via `availableStudies` and `availableGroups`
1. creates a connection to a specific study via `getStudy` or a group via `getGroup`
1. retrieves datasets by name via `getDataset`


### for example:

```{r}
library(DataSpaceR)

con <- connectDS()
con
```

`connectDS()` will create a connection to DataSpace.


### available studies can be listed by `availableStudies` field

```{r}
knitr::kable(head(con$availableStudies))
```


### available groups can be listed by `availableGroups` field

```{r}
knitr::kable(con$availableGroups)
```

***Note***: A group is a curated collection of participants from filtering of treatments, products, studies, or species, and it is created in [the DataSpace App](https://dataspace.cavd.org/cds/CAVD/app.view).

Check out [the reference page](https://docs.ropensci.org/DataSpaceR/reference/DataSpaceConnection.html) of `DataSpaceConnection` for all available fields and methods.


### create an instance of `cvd408`

```{r}
cvd408 <- con$getStudy("cvd408")
cvd408
class(cvd408)
```


### available datasets can be listed by `availableDatasets` field

```{r}
knitr::kable(cvd408$availableDatasets)
```

which will print names of available datasets.


### Neutralizing Antibody dataset (`NAb`) can be retrieved by:

```{r}
NAb <- cvd408$getDataset("NAb")
dim(NAb)
colnames(NAb)
```

Check out [the reference page](https://docs.ropensci.org/DataSpaceR/reference/DataSpaceStudy.html) of `DataSpaceStudy` for all available fields and methods.

***Note***: The package uses a [R6](https://cran.r-project.org/package=R6) class to represent the connection to a study and get around some of R's copy-on-change behavior.


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/DataSpaceR/issues).
* License: GPL-3
* Get citation information for `DataSpaceR` in R doing `citation(package = 'DataSpaceR')`
* Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/DataSpaceR/blob/main/CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Accessing Publication Data"
author: "Helen Miller"
date: "2021-08-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Accessing Publication Data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<style>
table {
    display: block;
    max-width: 100%;
    max-height: 600px;
    overflow: scroll;
}
thead{
    position: sticky;
    top: 0px;
    background-color: #fff;
}
</style>




DataSpace maintains a curated collection of relevant publications, which can be accessed through the [Publications page](https://dataspace.cavd.org/cds/CAVD/app.view?#learn/learn/Publication) through the app. Some publications laos include datasets which can be downloaded as a zip file. `DataSpaceR` provides an interface for browsing publications in DataSpace and downloading publication data where available.

## Browsing publications in DataSpace

The `DataSpaceConnection` object includes methods for browsing and downloading publications and publication data.


```r
library(DataSpaceR)
library(data.table)
con <- connectDS()
con
#> <DataSpaceConnection>
#>   URL: https://dataspace.cavd.org
#>   User: 
#>   Available studies: 273
#>     - 77 studies with data
#>     - 5049 subjects
#>     - 423195 data points
#>   Available groups: 3
#>   Available publications: 1530
#>     - 12 publications with data
```

The `DataSpaceConnection` print method summarizes the publications and publication data. More details about publications can be accessed through `con$availablePublications`.


```r
knitr::kable(head(con$availablePublications[, -"link"]))
```



|publication_id |first_author  |all_authors                                                                                                                                                                                                                                                              |title                                                                                                                                            |journal                     |publication_date |pubmed_id |related_studies |studies_with_data |publication_data_available |
|:--------------|:-------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------|:----------------|:---------|:---------------|:-----------------|:--------------------------|
|1006           |Abbink P      |Abbink P, Maxfield LF, Ng'ang'a D, Borducchi EN, Iampietro MJ, Bricault CA, Teigler JE, Blackmore S, Parenteau L, Wagh K, Handley SA, Zhao G, Virgin HW, Korber B, Barouch DH                                                                                            |Construction and evaluation of novel rhesus monkey adenovirus vaccine vectors                                                                    |J Virol                     |2015 Feb         |25410856  |NA              |NA                |FALSE                      |
|1303           |Abbott RK     |Abbott RK, Lee JH, Menis S, Skog P, Rossi M, Ota T, Kulp DW, Bhullar D, Kalyuzhniy O, Havenar-Daughton C, Schief WR, Nemazee D, Crotty S                                                                                                                                 |Precursor frequency and affinity determine B Cell competitive fitness in germinal centers, tested with germline-targeting HIV vaccine immunogens |Immunity                    |2018 Jan 16      |29287996  |NA              |NA                |FALSE                      |
|1136           |Abu-Raddad LJ |Abu-Raddad LJ, Boily MC, Self S, Longini IM Jr                                                                                                                                                                                                                           |Analytic insights into the population level impact of imperfect prophylactic HIV vaccines                                                        |J Acquir Immune Defic Syndr |2007 Aug 1       |17554215  |NA              |NA                |FALSE                      |
|1485           |Abu-Raddad LJ |Abu-Raddad LJ, Barnabas RV, Janes H, Weiss HA, Kublin JG, Longini IM Jr, Wasserheit JN, HIV Viral Load Working Group                                                                                                                                                     |Have the explosive HIV epidemics in sub-Saharan Africa been driven by higher community viral load?                                               |AIDS                        |2013 Mar 27      |23196933  |NA              |NA                |FALSE                      |
|842            |Acharya P     |Acharya P, Tolbert WD, Gohain N, Wu X, Yu L, Liu T, Huang W, Huang CC, Kwon YD, Louder RK, Luongo TS, McLellan JS, Pancera M, Yang Y, Zhang B, Flinko R, Foulke JS Jr, Sajadi MM, Kamin-Lewis R, Robinson JE, Martin L, Kwong PD, Guan Y, DeVico AL, Lewis GK, Pazgier M |Structural definition of an antibody-dependent cellular cytotoxicity response implicated in reduced risk for HIV-1 infection                     |J Virol                     |2014 Nov         |25165110  |NA              |NA                |FALSE                      |
|1067           |Ackerman ME   |Ackerman ME, Mikhailova A, Brown EP, Dowell KG, Walker BD, Bailey-Kellogg C, Suscovich TJ, Alter G                                                                                                                                                                       |Polyfunctional HIV-specific antibody responses are associated with spontaneous HIV control                                                       |PLOS Pathog                 |2016 Jan         |26745376  |NA              |NA                |FALSE                      |

This table summarizes all publications, providing some information like first author, journal where it was published, and title as a `data.table`. It also includes a pubmed url where available under `link`. Related studies under `related_studies`, and related studies with data available under `studies_with_data`.  We can use `data.table` methods to filter and sort this table to browse available publications.

For example, we can filter this table to view only publications related to a particular study:


```r
vtn096_pubs <- con$availablePublications[grepl("vtn096", related_studies)]
knitr::kable(vtn096_pubs[, -"link"])
```



|publication_id |first_author |all_authors                                                                                                                                                                                                                                                 |title                                                                                                                                                                              |journal              |publication_date |pubmed_id |related_studies        |studies_with_data |publication_data_available |
|:--------------|:------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------|:----------------|:---------|:----------------------|:-----------------|:--------------------------|
|1466           |Cram JA      |Cram JA, Fiore-Gartland AJ, Srinivasan S, Karuna S, Pantaleo G, Tomaras GD, Fredricks DN, Kublin JG                                                                                                                                                         |Human gut microbiota is associated with HIV-reactive immunoglobulin at baseline and following HIV vaccination                                                                      |PLoS One             |2019             |31869338  |vtn096                 |NA                |TRUE                       |
|250            |Huang Y      |Huang Y, DiazGranados C, Janes H, Huang Y, deCamp A, Metch B, Grant S, Sanchez B, Phogat S,  Koutsoukos M, Kanesa-Thasan N, Bourguignon P, Collard A, Buchbinder S, Tomaras GD, McElrath J, Gray G, Kublin J, Corey L, Gilbert PB                           |Selection of HIV vaccine candidates for concurrent testing in an efficacy trial                                                                                                    |Curr Opin Virol      |2016 Jan 28      |26827165  |mrv144, vtn096         |NA                |FALSE                      |
|267            |Huang Y      |Huang Y, Zhang L, Janes H, Frahm N, Isaacs A, Kim JH, Montefiori D, McElrath MJ, Tomaras GD, Gilbert PB                                                                                                                                                     |Predictors of durable immune responses six months after the last vaccination in preventive HIV vaccine trials                                                                      |Vaccine              |2017 Feb 22      |28131393  |mrv144, vtn096, vtn205 |NA                |FALSE                      |
|268            |Huang Y      |Huang Y, Gilbert P, Fu R, Janes H                                                                                                                                                                                                                           |Statistical methods for down-selection of treatment regimens based on multiple endpoints, with application to HIV vaccine trials                                                   |Biostatistics        |2017 Apr 1       |27649715  |mrv144, vtn096         |NA                |FALSE                      |
|1392           |Pantaleo G   |Pantaleo G, Janes H, Karuna S, Grant S, Ouedraogo GL, Allen M, Tomaras GD, Frahm N, Montefiori DC, Ferrari G, Ding S, Lee C, Robb ML, Esteban M, Wagner R, Bart PA, Rettby N, McElrath MJ, Gilbert PB, Kublin JG, Corey L, NIAID HIV Vaccine Trials Network |Safety and immunogenicity of a multivalent HIV vaccine comprising envelope protein with either DNA or NYVAC vectors (HVTN 096): a phase 1b, double-blind, placebo-controlled trial |Lancet HIV           |2019 Oct 7       |31601541  |mrv144, vtn096, vtn100 |NA                |TRUE                       |
|1420           |Westling T   |Westling T, Juraska M, Seaton KE, Tomaras GD, Gilbert PB, Janes H                                                                                                                                                                                           |Methods for comparing durability of immune responses between vaccine regimens in early-phase trials                                                                                |Stat Methods Med Res |2019 Jan 9       |30623732  |vtn094, vtn096         |NA                |FALSE                      |
|281            |Yates NL     |Yates NL, deCamp AC, Korber BT, Liao HX, Irene C, Pinter A, Peacock J, Harris LJ, Sawant S, Hraber P, Shen X, Rerks-Ngarm S, Pitisuttithum P, Nitayapan S, Berman PW, Robb ML, Pantaleo G, Zolla-Pazner S, Haynes BF, Alam SM, Montefiori DC, Tomaras GD    |HIV-1 envelope glycoproteins from diverse clades differentiate antibody responses and durability among vaccinees                                                                   |J Virol              |2018 Mar 28      |29386288  |mrv144, vtn096         |NA                |FALSE                      |


or publications that have related studies with integrated data in DataSpace:


```r
pubs_with_study_data <- con$availablePublications[!is.na(studies_with_data)]
knitr::kable(head(pubs_with_study_data[, -"link"]))
```



|publication_id |first_author      |all_authors                                                                                                                                                                                                                                                                                                                                                         |title                                                                                                                                         |journal            |publication_date |pubmed_id |related_studies                                                                |studies_with_data |publication_data_available |
|:--------------|:-----------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------|:------------------|:----------------|:---------|:------------------------------------------------------------------------------|:-----------------|:--------------------------|
|1540           |Andersen-Nissen E |Andersen-Nissen E, Fiore-Gartland A, Ballweber Fleming L, Carpp LN, Naidoo AF, Harper MS, Voillet V, Grunenberg N, Laher F, Innes C, Bekker L-G, Kublin JG, Huang Y, Ferrari G, Tomaras GD, Gray G, Gilbert PB, McElrath MJ                                                                                                                                         |Innate immune signatures to a partially-efficacious HIV vaccine predict correlates of HIV-1 infection risk                                    |PLOS Pathog        |2021             |NA        |vtn097                                                                         |vtn097            |TRUE                       |
|213            |Andrasik MP       |Andrasik MP, Yoon R, Mooney J, Broder G, Bolton M, Votto T, Davis-Vogel A; HVTN 505 study team; NIAID HIV Vaccine Trials Network                                                                                                                                                                                                                                    |Exploring barriers and facilitators to participation of male-to-female transgender persons in preventive HIV vaccine clinical trials.         |Prev Sci           |2014 Jun         |23446435  |vtn505                                                                         |vtn505            |FALSE                      |
|218            |Andrasik MP       |Andrasik MP, Karuna ST, Nebergall M, Koblin BA, Kublin JG; NIAID HIV Vaccine Trials Network                                                                                                                                                                                                                                                                         |Behavioral risk assessment in HIV Vaccine Trials Network (HVTN) clinical trials: A qualitative study exploring HVTN staff perspectives        |Vaccine            |2013 Sep 13      |23859840  |vtn069, vtn404, vtn502, vtn503, vtn504, vtn505, vtn802, vtn903, vtn906, vtn907 |vtn505            |FALSE                      |
|225            |Andrasik MP       |Andrasik MP, Chandler C, Powell B, Humes D, Wakefield S, Kripke K, Eckstein D                                                                                                                                                                                                                                                                                       |Bridging the divide: HIV prevention research and Black men who have sex with men                                                              |Am J Public Health |2014 Apr         |24524520  |vtn505                                                                         |vtn505            |FALSE                      |
|226            |Arnold MP         |Arnold MP, Andrasik M, Landers S, Karuna S, Mimiaga MJ, Wakefield S, Mayer K, Buchbinder S, Koblin BA                                                                                                                                                                                                                                                               |Sources of racial/ethnic differences in HIV vaccine trial awareness: Exposure, attention, or both?                                            |Am J Public Health |2014 Aug         |24922153  |vtn505                                                                         |vtn505            |FALSE                      |
|7              |Asbach B          |Asbach B, Kliche A, Köstler J, Perdiguero B, Esteban M, Jacobs BL, Montefiori DC, LaBranche CC, Yates NL, Tomaras GD, Ferrari G, Foulds KE, Roederer M, Landucci G, Forthal DN, Seaman MS, Hawkins N, Self SG, Sato A, Gottardo R, Phogat S, Tartaglia J, Barnett SW, Burke B, Cristillo AD, Weiss DE, Francis J, Galmin L, Ding S, Heeney JL, Pantaleo G, Wagner R |Potential to streamline heterologous DNA prime and NYVAC/protein boost HIV vaccine regimens in rhesus macaques by employing improved antigens |J Virol            |2016 Mar 28      |26865719  |cvd277, cvd281                                                                 |cvd277, cvd281    |FALSE                      |

We can also use this information to connect to related studies and pull integrated data. Say we are interested in this Rouphael (2019) publication:


```r
rouphael2019 <- con$availablePublications[first_author == "Rouphael NG"]
knitr::kable(rouphael2019[, -"link"])
```



|publication_id |first_author |all_authors                                                                                                                                                                                                                                                                                                                |title                                                                                         |journal       |publication_date |pubmed_id |related_studies |studies_with_data |publication_data_available |
|:--------------|:------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------|:-------------|:----------------|:---------|:---------------|:-----------------|:--------------------------|
|1390           |Rouphael NG  |Rouphael NG, Morgan C, Li SS, Jensen R, Sanchez B, Karuna S, Swann E, Sobieszczyk ME, Frank I, Wilson GJ, Tieu HV, Maenza J, Norwood A, Kobie J, Sinangil F, Pantaleo G, Ding S, McElrath MJ, De Rosa SC, Montefiori DC, Ferrari G, Tomaras GD, Keefer MC, HVTN 105 Protocol Team and the NIAID HIV Vaccine Trials Network |DNA priming and gp120 boosting induces HIV-specific antibodies in a randomized clinical trial |J Clin Invest |2019 Sep 30      |31566579  |vtn105          |vtn105            |TRUE                       |

We can find this publication in the available publications table, determine related studies, and pull data for those studies where available.


```r
related_studies <- rouphael2019$related_studies
related_studies
#> [1] "vtn105"
rouphael2019_study <- con$getStudy(related_studies)
rouphael2019_study
#> <DataSpaceStudy>
#>   Study: vtn105
#>   URL: https://dataspace.cavd.org/CAVD/vtn105
#>   Available datasets:
#>     - Binding Ab multiplex assay
#>     - Demographics
#>     - Intracellular Cytokine Staining
#>     - Neutralizing antibody
#>   Available non-integrated datasets:
dim(rouphael2019_study$availableDatasets)
#> [1] 4 4
```

We can see that there are datasets available for this study. We can pull any of them using `rouphael2019_study$getDataset()`.

## Downloading Publication Data


DataSpace also includes publication datasets for some publications. The format of this data will vary from publication to publication, and is stored in a zip file. The `publication_data_available` field specifies publications where publication data is available.


```r
pubs_with_data <- con$availablePublications[publication_data_available == TRUE]
knitr::kable(head(pubs_with_data[, -"link"]))
```



|publication_id |first_author      |all_authors                                                                                                                                                                                                                                                                                                                                                                                                                               |title                                                                                                                                                                              |journal         |publication_date |pubmed_id |related_studies |studies_with_data |publication_data_available |
|:--------------|:-----------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------|:----------------|:---------|:---------------|:-----------------|:--------------------------|
|1540           |Andersen-Nissen E |Andersen-Nissen E, Fiore-Gartland A, Ballweber Fleming L, Carpp LN, Naidoo AF, Harper MS, Voillet V, Grunenberg N, Laher F, Innes C, Bekker L-G, Kublin JG, Huang Y, Ferrari G, Tomaras GD, Gray G, Gilbert PB, McElrath MJ                                                                                                                                                                                                               |Innate immune signatures to a partially-efficacious HIV vaccine predict correlates of HIV-1 infection risk                                                                         |PLOS Pathog     |2021             |NA        |vtn097          |vtn097            |TRUE                       |
|136            |Bekker LG         |Bekker LG, Moodie Z, Grunenberg N, Laher F, Tomaras GD, Cohen KW, Allen M, Malahleha M, Mngadi K, Daniels B, Innes C, Bentley C, Frahm N, Morris DE, Morris L, Mkhize NN, Montefiori DC, Sarzotti-Kelsoe M, Grant S, Yu C, Mehra VL, Pensiero MN, Phogat S, DiazGranados CA, Barnett SW, Kanesa-thasan N, Koutsoukos M, Michael NL, Robb ML, Kublin JG, Gilbert PB, Corey L, Gray GE, McElrath MJ on behalf of the HVTN 100 Protocol Team |A phase 1/2 HIV-1 vaccine trial of a Subtype C ALVAC-HIV and Bivalent Subtype C gp120/MF59 vaccine regimen in low risk HIV uninfected South African adults                         |Lancet HIV      |2018 Jul         |29898870  |mrv144, vtn100  |NA                |TRUE                       |
|1466           |Cram JA           |Cram JA, Fiore-Gartland AJ, Srinivasan S, Karuna S, Pantaleo G, Tomaras GD, Fredricks DN, Kublin JG                                                                                                                                                                                                                                                                                                                                       |Human gut microbiota is associated with HIV-reactive immunoglobulin at baseline and following HIV vaccination                                                                      |PLoS One        |2019             |31869338  |vtn096          |NA                |TRUE                       |
|80             |Fong Y            |Fong Y, Shen X, Ashley VC, Deal A, Seaton KE, Yu C, Grant SP, Ferrari G, deCamp AC, Bailer RT, Koup RA, Montefiori D, Haynes BF, Sarzotti-Kelsoe M, Graham BS, Carpp LN, Hammer SM, Sobieszczyk M, Karuna S, Swann E, DeJesus E, Mulligan M, Frank I, Buchbinder S, Novak RM, McElrath MJ, Kalams S, Keefer M, Frahm NA, Janes HE, Gilbert PB, Tomaras GD                                                                                 |Modification of the association between T-Cell immune responses and human immunodeficiency virus type 1 infection risk by vaccine-induced antibody responses in the HVTN 505 trial |J Infect Dis    |2018 Mar 28      |29325070  |vtn505          |vtn505            |TRUE                       |
|79             |Hammer SM         |Hammer SC, Sobieszczyk ME, Janes H, Karuna ST, Mulligan MJ, Grove D, Koblin BA, Buchbinder SP, Keefer MC, Tomaras GD, Frahm N, Hural J, Anude C, Graham BS, Enama ME, Adams E, DeJesus E, Novak RM, Frank I, Bentley C, Ramirez S, Fu R, Koup RA, Mascola JR, Nabel GJ, Montefiori DC, Kublin J, McElrath MJ, Corey L, Gilbert PB                                                                                                         |Efficacy trial of a DNA/rAd5 HIV-1 preventive vaccine                                                                                                                              |N Engl J Med    |2013 Nov 28      |24099601  |vtn505          |vtn505            |TRUE                       |
|1467           |Hosseinipour MC   |Hosseinipour MC, Innes C, Naidoo S, Mann P, Hutter J, Ramjee G, Sebe M, Maganga L, Herce ME, deCamp AC, Marshall K, Dintwe O, Andersen-Nissen E, Tomaras GD, Mkhize N, Morris L, Jensen R, Miner MD, Pantaleo G, Ding S, Van Der Meeren O, Barnett SW, McElrath MJ, Corey L, Kublin JG                                                                                                                                                    |Phase 1 HIV vaccine trial to evaluate the safety and immunogenicity of HIV subtype C DNA and MF59-adjuvanted subtype C Env protein                                                 |Clin Infect Dis |2020 Jan 4       |31900486  |vtn111          |NA                |TRUE                       |

Data for a publication can be accessed through `DataSpaceR` with `con$downloadPublicationData()`. The publication ID must be specified, as found under `publication_id` in `con$availablePublications`. The file is presented as a zip file. The `unzip` argument gives us the option whether to unzip this file. By default, the file will be unzipped. You may also specify the directory to download the file. By default, it will be saved to your `Downloads` directory. This function invisibly returns the paths to the downloaded files.

Here, we download data for publication with ID 1461 (Westling, 2020), and view the resulting downloads.


```r
publicationDataFiles <- con$downloadPublicationData("1461", outputDir = tempdir(), unzip = TRUE, verbose = TRUE)
basename(publicationDataFiles)
#> [1] "causal.isoreg.fns.R"           "CD.SuperLearner.R"            
#> [3] "cd4_analysis.R"                "cd4_data.csv"                 
#> [5] "cd8_analysis.R"                "cd8_data.csv"                 
#> [7] "README.txt"                    "Westling_1461_file_format.pdf"
```

All zip files will include a file format document as a PDF, as well as a README. These documents will give an overview of the remaining contents of the files. In this case, data is separated by CD8+ T-cell responses and CD4+ T-cell responses, as described in the `README.txt`.


```r
cd4 <- fread(publicationDataFiles[grepl("cd4_data", publicationDataFiles)])
cd4
#>       pub_id age    sex      bmi   prot num_vacc dose response sexFemale
#>   1: 069-071  23   Male 33.31000 vtn069        3  4.0        0         0
#>   2: 069-068  20 Female 32.74000 vtn069        3  4.0        1         1
#>   3: 069-020  19   Male 20.24000 vtn069        3  4.0        0         0
#>   4: 069-008  19   Male 29.98000 vtn069        3  4.0        0         0
#>   5: 069-058  27   Male 31.48000 vtn069        3  4.0        0         0
#>  ---                                                                    
#> 368: 100-182  36 Female 24.32323 vtn100        4  1.5        0         1
#> 369: 100-034  20   Male 25.15590 vtn100        4  1.5        1         0
#> 370: 100-115  20   Male 22.94812 vtn100        4  1.5        1         0
#> 371: 100-144  19   Male 20.95661 vtn100        4  1.5        0         0
#> 372: 100-189  24   Male 19.03114 vtn100        4  1.5        1         0
#>      studyHVTN052 studyHVTN068 studyHVTN069 studyHVTN204 studyHVTN100 vacc_type
#>   1:            0            0            1            0            0   VRC4or6
#>   2:            0            0            1            0            0   VRC4or6
#>   3:            0            0            1            0            0   VRC4or6
#>   4:            0            0            1            0            0   VRC4or6
#>   5:            0            0            1            0            0   VRC4or6
#>  ---                                                                           
#> 368:            0            0            0            0            1   VRC4or6
#> 369:            0            0            0            0            1   VRC4or6
#> 370:            0            0            0            0            1   VRC4or6
#> 371:            0            0            0            0            1   VRC4or6
#> 372:            0            0            0            0            1   VRC4or6
#>      vacc_typeSinglePlasmid vacc_typeVRC4or6
#>   1:                      0                1
#>   2:                      0                1
#>   3:                      0                1
#>   4:                      0                1
#>   5:                      0                1
#>  ---                                        
#> 368:                      0                1
#> 369:                      0                1
#> 370:                      0                1
#> 371:                      0                1
#> 372:                      0                1
```

This publication also includes analysis scripts used for the publication, which can allow users to reproduce the analysis and results.


## Session information


```r
sessionInfo()
#> R version 4.1.0 (2021-05-18)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 20.04.2 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] data.table_1.14.0 DataSpaceR_0.7.5  knitr_1.33       
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.7       digest_0.6.27    assertthat_0.2.1 R6_2.5.0        
#>  [5] jsonlite_1.7.2   magrittr_2.0.1   evaluate_0.14    highr_0.9       
#>  [9] httr_1.4.2       stringi_1.7.3    curl_4.3.2       tools_4.1.0     
#> [13] stringr_1.4.0    Rlabkey_2.8.0    xfun_0.25        compiler_4.1.0
```
---
title: "Accessing Non-Integrated Datasets"
author: "Helen Miller"
date: "2021-08-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Accessing Non-Integrated Datasets}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



Many studies include data from assays which have not been integrated into the DataSpace. Some of these are available as "Non-Integrated Datasets," which can be downloaded from the app as a zip file. `DataSpaceR` provides an interface for accessing non-integrated data from studies where it is available.

## Viewing available non-integrated data

Methods on the DataSpace Study object allow you to see what non-integrated data may be available before downloading it. We will be using HVTN 505 as an example.


```r
library(DataSpaceR)
con <- connectDS()
vtn505 <- con$getStudy("vtn505")
vtn505
#> <DataSpaceStudy>
#>   Study: vtn505
#>   URL: https://dataspace.cavd.org/CAVD/vtn505
#>   Available datasets:
#>     - Binding Ab multiplex assay
#>     - Demographics
#>     - Intracellular Cytokine Staining
#>     - Neutralizing antibody
#>   Available non-integrated datasets:
#>     - ADCP
#>     - Demographics (Supplemental)
#>     - Fc Array
```

The print method on the study object will list available non-integrated datasets. The `availableDatasets` property shows some more info about available datasets, with the `integrated` field indicating whether the data is integrated. The value for `n` will be `NA` for non-integrated data until the dataset has been loaded.


```r
knitr::kable(vtn505$availableDatasets)
```



|name         |label                           |     n|integrated |
|:------------|:-------------------------------|-----:|:----------|
|BAMA         |Binding Ab multiplex assay      | 10260|TRUE       |
|Demographics |Demographics                    |  2504|TRUE       |
|ICS          |Intracellular Cytokine Staining | 22684|TRUE       |
|NAb          |Neutralizing antibody           |   628|TRUE       |
|ADCP         |ADCP                            |    NA|FALSE      |
|DEM SUPP     |Demographics (Supplemental)     |    NA|FALSE      |
|Fc Array     |Fc Array                        |    NA|FALSE      |

## Loading non-integrated data

Non-Integrated datasets can be loaded with `getDataset` like integrated data. This will unzip the non-integrated data to a temp directory and load it into the environment.


```r
adcp <- vtn505$getDataset("ADCP")
#> downloading vtn505_adcp.zip to /var/folders/d9/q394nzdj7b7_5hqkp1bvh9v40000gn/T//Rtmp4FHYQc...
#> No encoding supplied: defaulting to UTF-8.
#> unzipping vtn505_adcp.zip to /var/folders/d9/q394nzdj7b7_5hqkp1bvh9v40000gn/T//Rtmp4FHYQc/vtn505_adcp
dim(adcp)
#> [1] 378  11
colnames(adcp)
#>  [1] "study_prot"             "participant_id"         "study_day"              "lab_code"               "specimen_type"         
#>  [6] "antigen"                "percent_cv"             "avg_phagocytosis_score" "positivity_threshold"   "response"              
#> [11] "assay_identifier"
```

You can also view the file format info using `getDatasetDescription`. For non-integrated data, this will open a pdf into your computer's default pdf viewer.


```r
vtn505$getDatasetDescription("ADCP")
```

Non-integrated data is downloaded to a temp directory by default. There are a couple of ways to override this if desired. One is to specify `outputDir` when calling `getDataset` or `getDatasetDescription`.

If you will be accessing the data at another time and don't want to have to re-download it, you can change the default directory for the whole study object with `setDataDir`.


```r
vtn505$dataDir
#> [1] "/var/folders/d9/q394nzdj7b7_5hqkp1bvh9v40000gn/T//Rtmp4FHYQc"
vtn505$setDataDir(".")
vtn505$dataDir
#> [1] "/Users/juyeongkim/git/CDS/DataSpaceR/vignettes"
```

If the dataset already exists in the specified `dataDir` or `outputDir`, it will be not be downloaded. This can be overridden with `reload=TRUE`, which forces a re-download.


## Session information


```r
sessionInfo()
#> R version 4.1.1 (2021-08-10)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 11.5.2
#> 
#> Matrix products: default
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] DataSpaceR_0.7.5 knitr_1.33      
#> 
#> loaded via a namespace (and not attached):
#>  [1] httr_1.4.2        pkgload_1.2.1     jsonlite_1.7.2    spelling_2.2      assertthat_0.2.1  askpass_1.1       highr_0.9        
#>  [8] triebeard_0.3.0   urltools_1.7.3    yaml_2.2.1        remotes_2.4.0     sessioninfo_1.1.1 pillar_1.6.2      glue_1.4.2       
#> [15] uuid_0.1-4        digest_0.6.27     htmltools_0.5.2   pkgconfig_2.0.3   devtools_2.4.2    rcmdcheck_1.3.3   httpcode_0.3.0   
#> [22] purrr_0.3.4       gitcreds_0.1.1    codemetar_0.3.2   processx_3.5.2    tibble_3.1.4      openssl_1.4.4     usethis_2.0.1    
#> [29] ellipsis_0.3.2    cachem_1.0.6      withr_2.4.2       credentials_1.3.1 lazyeval_0.2.2    cli_3.0.1         magrittr_2.0.1   
#> [36] crayon_1.4.1      memoise_2.0.0     evaluate_0.14     ps_1.6.0          fs_1.5.0          fansi_0.5.0       xml2_1.3.2       
#> [43] parsedate_1.2.1   pkgbuild_1.2.0    tools_4.1.1       gh_1.3.0          hunspell_3.0.1    data.table_1.14.0 prettyunits_1.1.1
#> [50] Rlabkey_2.8.0     lifecycle_1.0.0   gert_1.3.2        stringr_1.4.0     xopen_1.0.0       pingr_2.0.1       callr_3.7.0      
#> [57] rex_1.2.0         compiler_4.1.1    covr_3.5.1        rhub_1.1.1        rlang_0.4.11      rstudioapi_0.13   sys_3.4          
#> [64] rappdirs_0.3.3    rmarkdown_2.10    testthat_3.0.4    curl_4.3.2        rematch_1.0.1     rematch2_2.1.2    R6_2.5.1         
#> [71] fastmap_1.1.0     utf8_1.2.2        commonmark_1.7    rprojroot_2.0.2   desc_1.3.0        stringi_1.7.4     crul_1.1.0       
#> [78] whoami_1.3.0      Rcpp_1.0.7        vctrs_0.3.8       xfun_0.25
```
---
title: "Accessing Monoclonal Antibody Data"
author:
- Ju Yeong Kim
- Jason Taylor
date: "2021-08-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Accessing Monoclonal Antibody Data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<style>
table {
    display: block;
    max-width: 100%;
    max-height: 600px;
    overflow: scroll;
}
thead{
    position: sticky;
    top: 0px;
    background-color: #fff;
}
</style>




## Workflow overview

Using the DataSpace [app](https://dataspace.cavd.org/cds/CAVD/app.view#mabgrid), the workflow of using the mAb grid is the following:

1. Navigate to the mAb Grid and browse the available mAb mixtures
2. Select the mAb mixtures that you'd like to investigate
3. Or filter rows by using columns:
    - mAb/Mixture
    - donor species
    - isotype
    - HXB2 location
    - tiers
    - clades
    - viruses
4. Click "Neutralization Curves" or "IC50 Titer Heatmap" to visualize the mAb data
5. Click "Export CSV" or "Export Excel" to download the mAb data

`DataSpaceR` offers a similar interface:

1. Browse the mAb Grid by `con$mabGridSummary`
2. Select the mAb mixtures by filtering the mAb grid using any columns found in `con$mabGrid` using `con$filterMabGrid()`
3. Use `con$getMab()` to retrieve the mAb data


## Browse the mAb Grid

You can browse the mAb Grid by calling the `mabGridSummary` field in the connection object:


```r
library(DataSpaceR)
con <- connectDS()

knitr::kable(head(con$mabGridSummary))
```



|mab_mixture    |donor_species |isotype |hxb2_location | n_viruses| n_clades| n_tiers| geometric_mean_curve_ic50| n_studies|
|:--------------|:-------------|:-------|:-------------|---------:|--------:|-------:|-------------------------:|---------:|
|10-1074        |human         |IgG     |Env           |         7|        3|       2|                 0.0213723|         1|
|10E8           |human         |IgG3    |gp160         |       227|       11|       7|                 0.4843333|         2|
|10E8 V2.0      |human         |        |              |        28|        3|       1|                 0.0031350|         1|
|10E8 V2.0/iMab |human         |        |gp160         |        13|        8|       2|                 0.0462897|         1|
|10E8 V4.0      |human         |        |              |        28|        3|       1|                 0.0024094|         1|
|10E8 V4.0/iMab |human         |        |              |       119|       12|       6|                 0.0015396|         1|

This table is designed to mimic the mAb grid found in the app.

One can also access the unsummarized data from the mAb grid by calling `con$mabGrid`.

## Filter the mAb grid

You can filter rows in the grid by specifying the values to keep in the columns found in the field `con$mabGrid`: `mab_mixture`, `donor_species`, `isotype`, `hxb2_location`, `tiers`, `clades`, `viruses`, and `studies`. `filterMabGrid` takes the column and the values and filters the underlying tables (private fields), and when you call the `mabGridSummary` or (which is actually an [active binding](https://r6.r-lib.org/articles/Introduction.html#active-bindings)), it returns the filtered grid with updated `n_` columns and `geometric_mean_curve_ic50`.


```r
# filter the grid by viruses
con$filterMabGrid(using = "virus", value = c("242-14", "Q23.17", "6535.3", "BaL.26", "DJ263.8"))

# filter the grid by donor species (llama)
con$filterMabGrid(using = "donor_species", value = "llama")

# check the updated grid
knitr::kable(con$mabGridSummary)
```



|mab_mixture |donor_species |isotype |hxb2_location | n_viruses| n_clades| n_tiers| geometric_mean_curve_ic50| n_studies|
|:-----------|:-------------|:-------|:-------------|---------:|--------:|-------:|-------------------------:|---------:|
|11F1B       |llama         |        |              |         4|        2|       1|                        NA|         1|
|11F1F       |llama         |        |              |         4|        2|       1|                26.2961178|         1|
|1H9         |llama         |        |Env           |         4|        2|       1|                 5.0898322|         1|
|2B4F        |llama         |        |              |         4|        2|       1|                 1.5242288|         1|
|2H10        |llama         |        |              |         2|        2|       1|                        NA|         1|
|2H10/W100A  |llama         |        |              |         2|        2|       1|                        NA|         1|
|3E3         |llama         |        |gp160         |         3|        3|       1|                 0.9944945|         1|
|4H73        |llama         |        |              |         4|        2|       1|                        NA|         1|
|5B10D       |llama         |        |              |         4|        2|       1|                        NA|         1|
|9B6B        |llama         |        |              |         4|        2|       1|                24.3643637|         1|
|A14         |llama         |        |gp160         |         3|        3|       1|                 1.8444582|         1|
|B21         |llama         |        |gp160         |         3|        3|       1|                 0.0936399|         1|
|B9          |llama         |        |gp160         |         3|        3|       1|                 0.0386986|         1|
|LAB5        |llama         |        |              |         4|        2|       1|                        NA|         1|

Or we can use method chaining to call multiple filter methods and browse the grid. Method chaining is unique to R6 objects and related to the pipe. See Hadley Wickham's [Advanced R](https://adv-r.hadley.nz/r6.html) for more info


```r
con$
  filterMabGrid(using = "hxb2_location", value = c("Env", "gp160"))$
  filterMabGrid(using = "donor_species", value = "llama")$
  mabGridSummary
```


## Retrieve column values from the mAb grid

You can retrieve values from the grid by `mab_mixture`, `donor_species`, `isotype`, `hxb2_location`, `tier`, `clade`, `virus`, and `studies`, or any variables found in the `mabGrid` field in the connection object via `data.table` operations.


```r
# retrieve available viruses in the filtered grid
con$mabGrid[, unique(virus)]
#> [1] "6535.3"  "Q23.17"  "242-14"  "BaL.26"  "DJ263.8"

# retrieve available clades for 1H9 mAb mixture in the filtered grid
con$mabGrid[mab_mixture == "1H9", unique(clade)]
#> [1] "B"        "CRF02_AG"
```


## Create a DataSpaceMab object

After filtering the grid, you can create a DataSpaceMab object that contains the filtered mAb data.


```r
mab <- con$getMab()
mab
#> <DataSpaceMab>
#>   URL: https://dataspace.cavd.org
#>   User: 
#>   Summary:
#>     - 3 studies
#>     - 14 mAb mixtures
#>     - 1 neutralization tiers
#>     - 3 clades
#>   Filters:
#>     - virus: 242-14, Q23.17, 6535.3, BaL.26, DJ263.8
#>     - mab_donor_species: llama
```

There are 6 public fields available in the `DataSpaceMab` object: `studyAndMabs`, `mabs`, `nabMab`, `studies`, `assays`, and `variableDefinitions`, and they are equivalent to the sheets in the excel file or the csv files you would download from the app via "Export Excel"/"Export CSV".


```r
knitr::kable(con$mabGridSummary)
```



|mab_mixture |donor_species |isotype |hxb2_location | n_viruses| n_clades| n_tiers| geometric_mean_curve_ic50| n_studies|
|:-----------|:-------------|:-------|:-------------|---------:|--------:|-------:|-------------------------:|---------:|
|11F1B       |llama         |        |              |         4|        2|       1|                        NA|         1|
|11F1F       |llama         |        |              |         4|        2|       1|                26.2961178|         1|
|1H9         |llama         |        |Env           |         4|        2|       1|                 5.0898322|         1|
|2B4F        |llama         |        |              |         4|        2|       1|                 1.5242288|         1|
|2H10        |llama         |        |              |         2|        2|       1|                        NA|         1|
|2H10/W100A  |llama         |        |              |         2|        2|       1|                        NA|         1|
|3E3         |llama         |        |gp160         |         3|        3|       1|                 0.9944945|         1|
|4H73        |llama         |        |              |         4|        2|       1|                        NA|         1|
|5B10D       |llama         |        |              |         4|        2|       1|                        NA|         1|
|9B6B        |llama         |        |              |         4|        2|       1|                24.3643637|         1|
|A14         |llama         |        |gp160         |         3|        3|       1|                 1.8444582|         1|
|B21         |llama         |        |gp160         |         3|        3|       1|                 0.0936399|         1|
|B9          |llama         |        |gp160         |         3|        3|       1|                 0.0386986|         1|
|LAB5        |llama         |        |              |         4|        2|       1|                        NA|         1|

## View metadata concerning the mAb object

There are several metadata fields that can be exported in the mAb object.


```r
names(mab)
#>  [1] ".__enclos_env__"     "variableDefinitions" "assays"             
#>  [4] "studies"             "nabMab"              "mabs"               
#>  [7] "studyAndMabs"        "config"              "clone"              
#> [10] "refresh"             "print"               "initialize"
```

## Session information


```r
sessionInfo()
#> R version 4.1.0 (2021-05-18)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 20.04.2 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] data.table_1.14.0 DataSpaceR_0.7.5  knitr_1.33       
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.7       digest_0.6.27    assertthat_0.2.1 R6_2.5.0        
#>  [5] jsonlite_1.7.2   magrittr_2.0.1   evaluate_0.14    highr_0.9       
#>  [9] httr_1.4.2       stringi_1.7.3    curl_4.3.2       tools_4.1.0     
#> [13] stringr_1.4.0    Rlabkey_2.8.0    xfun_0.25        compiler_4.1.0
```
---
title: "Introduction to DataSpaceR"
author: "Ju Yeong Kim"
date: "2021-08-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to DataSpaceR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<style>
table {
    display: block;
    max-width: 100%;
    max-height: 600px;
    overflow: scroll;
}
thead{
    position: sticky;
    top: 0px;
    background-color: #fff;
}
</style>



This package provides a thin wrapper around [Rlabkey](https://cran.r-project.org/package=Rlabkey) and connects to the the [CAVD DataSpace](https://dataspace.cavd.org) database, making it easier to fetch datasets from specific studies.

## Configuration

First, go to [DataSpace](https://dataspace.cavd.org) now and set yourself up with an account.

In order to connect to the CAVD DataSpace via `DataSpaceR`, you will need a `netrc` file in your home directory that will contain a `machine` name (hostname of DataSpace), and `login` and `password`. There are two ways to create a `netrc` file.

### Creating a netrc file with `writeNetrc`

On your R console, create a `netrc` file using a function from `DataSpaceR`:


```r
writeNetrc(
  login = "yourEmail@address.com",
  password = "yourSecretPassword",
  netrcFile = "/your/home/directory/.netrc" # use getNetrcPath() to get the default path
)
```

This will create a `netrc` file in your home directory. Make sure you have a valid login and password.

### Manually creating a netrc file

***Alternatively***, you can manually create a netrc file.

* On Windows, this file should be named `_netrc`
* On UNIX/Mac, it should be named `.netrc`
* The file should be located in the user's home directory, and the permissions on the file should be unreadable for everybody except the owner
* To determine your home directory, run `Sys.getenv("HOME")` in R

The following three lines must be included in the `.netrc` or `_netrc` file either separated by white space (spaces, tabs, or newlines) or commas. Multiple such blocks can exist in one file.

```
machine dataspace.cavd.org
login myuser@domain.com
password supersecretpassword
```

See [here](https://www.labkey.org/wiki/home/Documentation/page.view?name=netrc) for more information about `netrc`.


## Initiate a connection

We'll be looking at study `cvd256`. If you want to use a different study, change that string. You can instantiate multiple connections to different studies simultaneously.


```r
library(DataSpaceR)
con <- connectDS()
con
#> <DataSpaceConnection>
#>   URL: https://dataspace.cavd.org
#>   User: 
#>   Available studies: 273
#>     - 77 studies with data
#>     - 5049 subjects
#>     - 423195 data points
#>   Available groups: 3
#>   Available publications: 1530
#>     - 12 publications with data
```

The call to `connectDS` instantiates the connection. Printing the object shows where it's connected and the available studies.


```r
knitr::kable(head(con$availableStudies))
```



|study_name |short_name                   |title                                                                                                                                                                                                                                                               |type               |status   |stage            |species            |start_date |strategy                             |network |data_availability |ni_data_availability |
|:----------|:----------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:------------------|:--------|:----------------|:------------------|:----------|:------------------------------------|:-------|:-----------------|:--------------------|
|cor01      |NA                           |The correlate of risk targeted intervention study (CORTIS):  A randomized, partially-blinded, clinical trial of isoniazid and rifapentine (3HP) therapy to prevent pulmonary tuberculosis in high-risk individuals identified by a transcriptomic correlate of risk |Phase III          |Inactive |Assays Completed |Human              |NA         |NA                                   |GH-VAP  |NA                |NA                   |
|cvd232     |Parks_RV_232                 |​Limiting Dose Vaginal SIVmac239 Challenge of RhCMV-SIV vaccinated Indian rhesus macaques.                                                                                                                                                                           |Pre-Clinical NHP   |Inactive |Assays Completed |Rhesus macaque     |2009-11-24 |Vector vaccines (viral or bacterial) |CAVD    |NA                |NA                   |
|cvd234     |Zolla-Pazner_Mab_test1 Study |Zolla-Pazner_Mab_Test1                                                                                                                                                                                                                                              |Antibody Screening |Inactive |Assays Completed |Non-Organism Study |2009-02-03 |Prophylactic neutralizing Ab         |CAVD    |NA                |NA                   |
|cvd235     |mAbs potency                 |Weiss mAbs potency                                                                                                                                                                                                                                                  |Antibody Screening |Inactive |Assays Completed |Non-Organism Study |2008-08-21 |Prophylactic neutralizing Ab         |CAVD    |NA                |NA                   |
|cvd236     |neutralization assays        |neutralization assays                                                                                                                                                                                                                                               |Antibody Screening |Active   |In Progress      |Non-Organism Study |2009-02-03 |Prophylactic neutralizing Ab         |CAVD    |NA                |NA                   |
|cvd238     |Gallo_PA_238                 |HIV-1 neutralization responses in chronically infected individuals                                                                                                                                                                                                  |Antibody Screening |Inactive |Assays Completed |Non-Organism Study |2009-01-08 |Prophylactic neutralizing Ab         |CAVD    |NA                |NA                   |

`con$availableStudies` shows the available studies in the CAVD DataSpace. Check out [the reference page](https://docs.ropensci.org/DataSpaceR/reference/DataSpaceConnection.html) of `DataSpaceConnection` for all available fields and methods.


```r
cvd256 <- con$getStudy("cvd256")
cvd256
#> <DataSpaceStudy>
#>   Study: cvd256
#>   URL: https://dataspace.cavd.org/CAVD/cvd256
#>   Available datasets:
#>     - Binding Ab multiplex assay
#>     - Demographics
#>     - Neutralizing antibody
#>   Available non-integrated datasets:
```

`con$getStudy` creates a connection to the study `cvd256`. Printing the object shows where it's connected, to what study, and the available datasets.


```r
knitr::kable(cvd256$availableDatasets)
```



|name         |label                      |    n|integrated |
|:------------|:--------------------------|----:|:----------|
|BAMA         |Binding Ab multiplex assay | 6740|TRUE       |
|Demographics |Demographics               |  121|TRUE       |
|NAb          |Neutralizing antibody      | 1419|TRUE       |

```r
knitr::kable(cvd256$treatmentArm)
```



|arm_id        |arm_part |arm_group |arm_name |randomization |coded_label     | last_day|description                                                                                           |
|:-------------|:--------|:---------|:--------|:-------------|:---------------|--------:|:-----------------------------------------------------------------------------------------------------|
|cvd256-NA-A-A |NA       |A         |A        |Vaccine       |Group A Vaccine |      168|DNA-C 4 mg administered IM at weeks 0, 4, and 8 AND NYVAC-C 10^7pfu/mL administered IM at week 24     |
|cvd256-NA-B-B |NA       |B         |B        |Vaccine       |Group B Vaccine |      168|DNA-C 4 mg administered IM at weeks 0 and 4 AND NYVAC-C 10^7pfu/mL administered IM at weeks 20 and 24 |

Available datasets and treatment arm information for the connection can be accessed by `availableDatasets` and `treatmentArm`.


## Fetching datasets

We can grab any of the datasets listed in the connection (`availableDatasets`).


```r
NAb <- cvd256$getDataset("NAb")
dim(NAb)
#> [1] 1419   33
colnames(NAb)
#>  [1] "participant_id"      "participant_visit"   "visit_day"          
#>  [4] "assay_identifier"    "summary_level"       "specimen_type"      
#>  [7] "antigen"             "antigen_type"        "virus"              
#> [10] "virus_type"          "virus_insert_name"   "clade"              
#> [13] "neutralization_tier" "tier_clade_virus"    "target_cell"        
#> [16] "initial_dilution"    "titer_ic50"          "titer_ic80"         
#> [19] "response_call"       "nab_lab_source_key"  "lab_code"           
#> [22] "exp_assayid"         "titer_id50"          "titer_id80"         
#> [25] "nab_response_id50"   "nab_response_id80"   "slope"              
#> [28] "vaccine_matched"     "study_prot"          "virus_full_name"    
#> [31] "virus_species"       "virus_host_cell"     "virus_backbone"
```

The *cvd256* object is an [`R6`](https://cran.r-project.org/package=R6) class, so it behaves like a true object. Functions (like `getDataset`) are members of the object, thus the `$` semantics to access member functions.

We can get detailed variable information using `getDatasetDescription`. `getDataset` and `getDatasetDescription` accept either the `name` or `label` field listed in `availableDatasets`.


```r
knitr::kable(cvd256$getDatasetDescription("NAb"))
```



|fieldName           |caption                                     |type                 |description                                                                                                                                                               |
|:-------------------|:-------------------------------------------|:--------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|ParticipantId       |Participant ID                              |Text (String)        |Subject identifier                                                                                                                                                        |
|antigen             |Antigen name                                |Text (String)        |The name of the antigen (virus) being tested.                                                                                                                             |
|antigen_type        |Antigen type                                |Text (String)        |The standardized term for the type of virus used in the construction of the nAb antigen.                                                                                  |
|assay_identifier    |Assay identifier                            |Text (String)        |Name identifying assay                                                                                                                                                    |
|clade               |Virus clade                                 |Text (String)        |The clade (gene subtype) of the virus (antigen) being tested.                                                                                                             |
|exp_assayid         |Experimental Assay Design Code              |Integer              |Unique ID assigned to the experiment design of the assay for tracking purposes.                                                                                           |
|initial_dilution    |Initial dilution                            |Number (Double)      |Indicates the initial specimen dilution.                                                                                                                                  |
|lab_code            |Lab ID                                      |Text (String)        |A code indicating the lab performing the assay.                                                                                                                           |
|nab_lab_source_key  |Data provenance                             |Integer              |Details regarding the provenance of the assay results.                                                                                                                    |
|nab_response_ID50   |Response call ID50                          |True/False (Boolean) |Indicates if neutralization is detected based on ID50 titer.                                                                                                              |
|nab_response_ID80   |Response call ID80                          |True/False (Boolean) |Indicates if neutralization is detected based on ID80 titer.                                                                                                              |
|neutralization_tier |Neutralization tier                         |Text (String)        |A classification specific to HIV NAb assay design, in which an antigen is assessed for its ease of neutralization (1=most easily neutralized, 3=least easily neutralized) |
|response_call       |Response call                               |True/False (Boolean) |Indicates if neutralization is detected.                                                                                                                                  |
|slope               |Slope                                       |Number (Double)      |The slope calculated using the difference between 50% and 80% neutralization.                                                                                             |
|specimen_type       |Specimen type                               |Text (String)        |The type of specimen used in the assay. For nAb assays, this is generally serum or plasma.                                                                                |
|study_prot          |Study Protocol                              |Text (String)        |Study protocol                                                                                                                                                            |
|summary_level       |Data summary level                          |Text (String)        |Defines the level at which the magnitude or response has been summarized (e.g. summarized at the isolate level).                                                          |
|target_cell         |Target cell                                 |Text (String)        |The cell line used in the assay to determine infection (lack of neutralization).  Generally TZM-bl or A3R5, but can also be other cell lines or non-engineered cells.     |
|tier_clade_virus    |Neutralization tier + Antigen clade + Virus |Text (String)        |A combination of neutralization tier, antigen clade, and virus used for filtering.                                                                                        |
|titer_ID50          |Titer ID50                                  |Number (Double)      |The adjusted value of 50% maximal inhibitory dilution (ID50).                                                                                                             |
|titer_ID80          |Titer ID80                                  |Number (Double)      |The adjusted value of 80% maximal inhibitory dilution (ID80).                                                                                                             |
|titer_ic50          |Titer IC50                                  |Number (Double)      |The half maximal inhibitory concentration (IC50).                                                                                                                         |
|titer_ic80          |Titer IC80                                  |Number (Double)      |The 80% maximal inhibitory concentration (IC80).                                                                                                                          |
|vaccine_matched     |Antigen vaccine match indicator             |True/False (Boolean) |Indicates if the interactive part of the antigen was designed to match the immunogen in the vaccine.                                                                      |
|virus               |Virus name                                  |Text (String)        |The term for the virus (antigen) being tested.                                                                                                                            |
|virus_backbone      |Virus backbone                              |Text (String)        |Indicates the backbone used to generate the virus if from a different plasmid than the envelope.                                                                          |
|virus_full_name     |Virus full name                             |Text (String)        |The full name of the virus used in the construction of the nAb antigen.                                                                                                   |
|virus_host_cell     |Virus host cell                             |Text (String)        |The host cell used to incubate the virus stock.                                                                                                                           |
|virus_insert_name   |Virus insert name                           |Text (String)        |The amino acid sequence inserted in the virus construct.                                                                                                                  |
|virus_species       |Virus species                               |Text (String)        |A classification for virus species using informal taxonomy.                                                                                                               |
|virus_type          |Virus type                                  |Text (String)        |The type of virus used in the construction of the nAb antigen.                                                                                                            |
|visit_day           |Visit Day                                   |Integer              |Target study day defined for a study visit. Study days are relative to Day 0, where Day 0 is typically defined as enrollment and/or first injection.                      |

To get only a subset of the data and speed up the download, filters can be passed to `getDataset`. The filters are created using the `makeFilter` function of the `Rlabkey` package.


```r
cvd256Filter <- makeFilter(c("visit_day", "EQUAL", "0"))
NAb_day0 <- cvd256$getDataset("NAb", colFilter = cvd256Filter)
dim(NAb_day0)
#> [1] 709  33
```

See `?makeFilter` for more information on the syntax.


## Creating a connection to all studies

To fetch data from multiple studies, create a connection at the project level.


```r
cavd <- con$getStudy("")
```

This will instantiate a connection at the `CAVD` level. Most functions work cross study connections just like they do on single studies.

You can get a list of datasets available across all studies.


```r
cavd
#> <DataSpaceStudy>
#>   Study: CAVD
#>   URL: https://dataspace.cavd.org/CAVD
#>   Available datasets:
#>     - Binding Ab multiplex assay
#>     - Demographics
#>     - Enzyme-Linked ImmunoSpot
#>     - Intracellular Cytokine Staining
#>     - Neutralizing antibody
#>     - PK MAb
#>   Available non-integrated datasets:
knitr::kable(cavd$availableDatasets)
```



|name         |label                           |      n|integrated |
|:------------|:-------------------------------|------:|:----------|
|BAMA         |Binding Ab multiplex assay      | 170320|TRUE       |
|Demographics |Demographics                    |   5049|TRUE       |
|ELISPOT      |Enzyme-Linked ImmunoSpot        |   5610|TRUE       |
|ICS          |Intracellular Cytokine Staining | 195883|TRUE       |
|NAb          |Neutralizing antibody           |  51382|TRUE       |
|PKMAb        |PK MAb                          |   3217|TRUE       |

In all-study connection, `getDataset` will combine the requested datasets. Note that in most cases, the datasets will have too many subjects for quick data transfer, making filtering of the data a necessity. The `colFilter` argument can be used here, as described in the `getDataset` section.


```r
conFilter <- makeFilter(c("species", "EQUAL", "Human"))
human <- cavd$getDataset("Demographics", colFilter = conFilter)
dim(human)
#> [1] 3142   36
colnames(human)
#>  [1] "subject_id"                      "subject_visit"                  
#>  [3] "species"                         "subspecies"                     
#>  [5] "sexatbirth"                      "race"                           
#>  [7] "ethnicity"                       "country_enrollment"             
#>  [9] "circumcised_enrollment"          "bmi_enrollment"                 
#> [11] "agegroup_range"                  "agegroup_enrollment"            
#> [13] "age_enrollment"                  "study_label"                    
#> [15] "study_start_date"                "study_first_enr_date"           
#> [17] "study_fu_complete_date"          "study_public_date"              
#> [19] "study_network"                   "study_last_vaccination_day"     
#> [21] "study_type"                      "study_part"                     
#> [23] "study_group"                     "study_arm"                      
#> [25] "study_arm_summary"               "study_arm_coded_label"          
#> [27] "study_randomization"             "study_product_class_combination"
#> [29] "study_product_combination"       "study_short_name"               
#> [31] "study_grant_pi_name"             "study_strategy"                 
#> [33] "study_prot"                      "genderidentity"                 
#> [35] "studycohort"                     "bmi_category"
```

Check out [the reference page](https://docs.ropensci.org/DataSpaceR/reference/DataSpaceStudy.html) of `DataSpaceStudy` for all available fields and methods.


## Connect to a saved group

A group is a curated collection of participants from filtering of treatments, products, studies, or species, and it is created in [the DataSpace App](https://dataspace.cavd.org/cds/CAVD/app.view).

Let's say you are using the App to filter and visualize data and want to save them for later or explore in R with `DataSpaceR`. You can save a group by clicking the Save button on the Active Filter Panel.

We can browse available the saved groups or the curated groups by DataSpace Team via `availableGroups`.


```r
knitr::kable(con$availableGroups)
```



| group_id|label                              |original_label                     |description                                                                                                               |created_by |shared |   n|studies                        |
|--------:|:----------------------------------|:----------------------------------|:-------------------------------------------------------------------------------------------------------------------------|:----------|:------|---:|:------------------------------|
|      220|NYVAC durability comparison        |NYVAC_durability                   |Compare durability in 4 NHP studies using NYVAC-C (vP2010)  and NYVAC-KC-gp140 (ZM96) products.                           |ehenrich   |TRUE   |  78|cvd281, cvd434, cvd259, cvd277 |
|      228|HVTN 505 case control subjects     |HVTN 505 case control subjects     |Participants from HVTN 505 included in the case-control analysis                                                          |drienna    |TRUE   | 189|vtn505                         |
|      230|HVTN 505 polyfunctionality vs BAMA |HVTN 505 polyfunctionality vs BAMA |Compares ICS polyfunctionality (CD8+, Any Env) to BAMA mfi-delta (single Env antigen) in the HVTN 505 case control cohort |drienna    |TRUE   | 170|vtn505                         |

To fetch data from a saved group, create a connection at the project level with a group ID. For example, we can connect to the "NYVAC durability comparison" group which has group ID 220 by `getGroup`.


```r
nyvac <- con$getGroup(220)
nyvac
#> <DataSpaceStudy>
#>   Group: NYVAC durability comparison
#>   URL: https://dataspace.cavd.org/CAVD
#>   Available datasets:
#>     - Binding Ab multiplex assay
#>     - Demographics
#>     - Enzyme-Linked ImmunoSpot
#>     - Intracellular Cytokine Staining
#>     - Neutralizing antibody
#>   Available non-integrated datasets:
```

Retrieving a dataset is the same as before.


```r
NAb_nyvac <- nyvac$getDataset("NAb")
dim(NAb_nyvac)
#> [1] 4281   33
```

## Access Virus Metadata

DataSpace maintains metadata about all viruses used in Neutralizing Antibody (NAb) assays. This data can be accessed through the app on the [NAb antigen page](https://dataspace.cavd.org/cds/CAVD/app.view#learn/learn/Assay/NAB/antigens) and [NAb MAb antigen page](https://dataspace.cavd.org/cds/CAVD/app.view#learn/learn/Assay/NAB%20MAB/antigens).

We can access this metadata in `DataSpaceR` with `con$virusMetadata`:


```r
knitr::kable(head(con$virusMetadata))
```



|assay_identifier |virus        |virus_type     |neutralization_tier |clade |antigen_control |virus_full_name                |virus_name_other |virus_species |virus_host_cell |virus_backbone |panel_names          |
|:----------------|:------------|:--------------|:-------------------|:-----|:---------------|:------------------------------|:----------------|:-------------|:---------------|:--------------|:--------------------|
|NAB MAB          |0013095-2.11 |Env Pseudotype |2                   |NA    |0               |0013095-2.11 [SG3Δenv] 293T/17 |NA               |HIV           |293T/17         |SG3Δenv        |Tiered diverse panel |
|NAB MAB          |001428-2.42  |Env Pseudotype |2                   |C     |0               |001428-2.42 [SG3Δenv] 293T/17  |NA               |HIV           |293T/17         |SG3Δenv        |Tiered diverse panel |
|NAB MAB          |0041.v3.c18  |Env Pseudotype |2                   |C     |0               |0041.v3.c18 [SG3Δenv] 293T/17  |0041.V3.C18      |HIV           |293T/17         |SG3Δenv        |NA                   |
|NAB MAB          |0077.v1.c16  |Env Pseudotype |2                   |C     |0               |0077.v1.c16 [SG3Δenv] 293T/17  |0077.v1.c16      |HIV           |293T/17         |SG3Δenv        |NA                   |
|NAB              |00836-2.5    |Env Pseudotype |1B                  |C     |0               |00836-2.5 [SG3Δenv] 293T/17    |NA               |HIV           |293T/17         |SG3Δenv        |Tiered diverse panel |
|NAB MAB          |0260.v5.c1   |Env Pseudotype |2                   |A     |0               |0260.v5.c1 [SG3Δenv] 293T/17   |0260.V5.C1       |HIV           |293T/17         |SG3Δenv        |Tiered diverse panel |

## Access monoclonal antibody data

See other vignette for a tutorial on accessing monoclonal antibody data with `DataSpaceR`:


```r
vignette("Monoconal_Antibody_Data")
```

## Browse and Download Publication Data

DataSpace maintains a curated collection of relevant publications, which can be accessed through the [Publications page](https://dataspace.cavd.org/cds/CAVD/app.view?#learn/learn/Publication) through the app. Metadata about these publications can be accessed through `DataSpaceR` with `con$availablePublications`.

See Publication Data vignette for a tutorial on accessing publication data through DataSpaceR.

```r
vignette("Publication_Data")
```

## Reference Tables

The followings are the tables of all fields and methods that work on `DataSpaceConnection` and `DataSpaceStudy` objects and could be used as a quick reference.

### `DataSpaceConnection`

| Name | Description |
| --- | --- |
| `availableStudies` | The table of available studies. |
| `availableGroups` | The table of available groups. |
| `availablePublications` | The table of available publications. |
| `mabGrid` | The filtered mAb grid. |
| `mabGridSummary` | The summarized mAb grid with updated `n_` columns and `geometric_mean_curve_ic50`. |
| `virusMetadata` | Metadata about all viruses in the DataSpace. |
| `filterMabGrid` | Filter rows in the mAb grid by specifying the values to keep in the columns found in the `mabGrid` field. |
| `resetMabGrid` | Reset the mAb grid to the unfiltered state. |
| `getMab` | Create a `DataSpaceMab` object by filtered `mabGrid`. |
| `getStudy` | Create a `DataSpaceStudy` object by study. |
| `getGroup` | Create a `DataSpaceStudy` object by group. |
| `downloadPublicationData` | Download data from a chosen publication. |


### `DataSpaceStudy`

| Name | Description |
| --- | --- |
| `study` | The study name. |
| `group` | The group name. |
| `availableDatasets` | The table of datasets available in the study object. |
| `treatmentArm` | The table of treatment arm information for the connected study. Not available for all study connection. |
| `dataDir` | The default target directory for downloading non-integrated datasets. |
| `studyInfo` | Stores the information about the study. |
| `getDataset` | Get a dataset from the connection. |
| `getDatasetDescription` | Get variable information. |
| `setDataDir` | Set default target directory for downloading non-integrated datasets. |

### `DataSpaceMab`

| Name | Description |
| --- | --- |
| `studyAndMabs` | The table of available mAbs by study. |
| `mabs` | The table of available mAbs and their attributes. |
| `nabMab` | The table of mAbs and their neutralizing measurements against viruses. |
| `studies` | The table of available studies. |
| `assays` | The table of assay status by study. |
| `variableDefinitions` | The table of variable definitions. |


## Session information


```r
sessionInfo()
#> R version 4.1.0 (2021-05-18)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 20.04.2 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.9.0
#> LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.9.0
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] data.table_1.14.0 DataSpaceR_0.7.5  knitr_1.33       
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.7       digest_0.6.27    assertthat_0.2.1 R6_2.5.0        
#>  [5] jsonlite_1.7.2   magrittr_2.0.1   evaluate_0.14    highr_0.9       
#>  [9] httr_1.4.2       stringi_1.7.3    curl_4.3.2       tools_4.1.0     
#> [13] stringr_1.4.0    Rlabkey_2.8.0    xfun_0.25        compiler_4.1.0
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netrc.R
\name{checkNetrc}
\alias{checkNetrc}
\title{Check netrc file}
\usage{
checkNetrc(netrcFile = getNetrcPath(), onStaging = FALSE, verbose = TRUE)
}
\arguments{
\item{netrcFile}{A character. File path to netrc file to check.}

\item{onStaging}{A logical. Whether to check the staging server instead
of the production server.}

\item{verbose}{A logical. Whether to print the extra details for
troubleshooting.}
}
\value{
The name of the netrc file
}
\description{
Check that there is a netrc file with a valid entry for the
CAVD DataSpace.
}
\examples{
try(checkNetrc())
}
\seealso{
\code{\link{connectDS}} \code{\link{writeNetrc}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/connectDS.R
\name{connectDS}
\alias{connectDS}
\title{Create a connection to DataSpace}
\usage{
connectDS(login = NULL, password = NULL, verbose = FALSE, onStaging = FALSE)
}
\arguments{
\item{login}{A character. Optional argument. If there is no netrc
file a temporary one can be written by passing login and password of an
active DataSpace account.}

\item{password}{A character. Optional. The password for the selected
login.}

\item{verbose}{A logical. Whether to print the extra details for
troubleshooting.}

\item{onStaging}{A logical. Whether to connect to the staging server instead
of the production server.}
}
\value{
an instance of \code{DataSpaceConnection}
}
\description{
Constructor for \code{\link{DataSpaceConnection}}
}
\details{
Instantiates an \code{DataSpaceConnection}.
The constructor will try to take the values of the various \code{labkey.*}
parameters from the global environment. If they don't exist, it will use
default values. These are assigned to `options`, which are then used by the
\code{DataSpaceConnection} class.
}
\examples{
\dontrun{
con <- connectDS()
}

con <- try(connectDS())
if (inherits(con, "try-error")) {
  warning("Read README for more information on how to set up a .netrc file.")
}
}
\seealso{
\code{\link{DataSpaceConnection}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataSpaceMab.R
\name{DataSpaceMab}
\alias{DataSpaceMab}
\title{The DataSpaceMab class}
\description{
The DataSpaceMab class

The DataSpaceMab class
}
\section{Constructor}{

\code{DataSpaceConnection$getMab()}
}

\examples{
\dontrun{
# Create a connection (Initiate a DataSpaceConnection object)
con <- connectDS()

# Browse the mAb Grid
con$mabGridSummary

# Filter the grid by viruses
con$filterMabGrid(using = "virus", value = c("242-14", "Q23.17", "6535.3", "BaL.26", "DJ263.8"))

# Filter the grid by donor species (llama)
con$filterMabGrid(using = "donor_species", value = "llama")

# Check the updated grid
con$mabGridSummary

# Retrieve available viruses in the filtered grid
con$mabGrid[, unique(virus)]

# Retrieve available clades for 1H9 mAb mixture in the filtered grid
con$mabGrid[mab_mixture == "1H9", unique(clade)]

# Create a DataSpaceMab object that contains the filtered mAb data
mab <- con$getMab()
mab

# Inspect the `nabMab` field
mab$nabMab
}

}
\seealso{
\code{\link{connectDS}} \code{\link{DataSpaceConnection}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{config}}{A list. Stores configuration of the connection object such
as URL, path and username.}

\item{\code{studyAndMabs}}{A data.table. The table of available mAbs by study.}

\item{\code{mabs}}{A data.table. The table of available mAbs and their
attributes.}

\item{\code{nabMab}}{A data.table. The table of mAbs and their neutralizing
measurements against viruses.}

\item{\code{studies}}{A data.table. The table of available studies.}

\item{\code{assays}}{A data.table. The table of assay status by study.}

\item{\code{variableDefinitions}}{A data.table. The table of variable
definitions.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{DataSpaceMab$new()}}
\item \href{#method-print}{\code{DataSpaceMab$print()}}
\item \href{#method-refresh}{\code{DataSpaceMab$refresh()}}
\item \href{#method-clone}{\code{DataSpaceMab$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize \code{DataSpaceMab} object.
See \code{\link{DataSpaceConnection}}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceMab$new(mabMixture, filters, config)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{mabMixture}}{A character vector.}

\item{\code{filters}}{A list.}

\item{\code{config}}{A list.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
Print the \code{DataSpaceMab} object summary.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceMab$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-refresh"></a>}}
\if{latex}{\out{\hypertarget{method-refresh}{}}}
\subsection{Method \code{refresh()}}{
Refresh the \code{DataSpaceMab} object to update datasets.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceMab$refresh()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceMab$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/DataSpaceR.R
\name{DataSpaceR-package}
\alias{DataSpaceR-package}
\alias{DataSpaceR}
\title{DataSpaceR}
\description{
DataSpaceR provides a convenient API for accessing datasets
within the DataSpace database.
}
\details{
Uses the Rlabkey package to connect to DataSpace. Implements
convenient methods for accessing datasets.
}
\seealso{
\code{\link{connectDS}}
}
\author{
Ju Yeong Kim
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataSpaceStudy.R
\name{DataSpaceStudy}
\alias{DataSpaceStudy}
\title{The DataSpaceStudy class}
\description{
The DataSpaceStudy class

The DataSpaceStudy class
}
\section{Constructor}{

\code{DataSpaceConnection$getStudy()}
\code{DataSpaceConnection$getGroup()}
}

\examples{
\dontrun{
# Create a connection (Initiate a DataSpaceConnection object)
con <- connectDS()

# Connect to cvd408 (Initiate a DataSpaceStudy object)
# https://dataspace.cavd.org/cds/CAVD/app.view#learn/learn/Study/cvd408?q=408
cvd408 <- con$getStudy("cvd408")
cvd408

# Retrieve Neutralizing antibody dataset (NAb) for cvd408 from DataSpace
NAb <- cvd408$getDataset("NAb")

# Get variable information of the NAb dataset
cvd408$getDatasetDescription("NAb")

# Take a look at cvd408's treatment arm information
cvd408$treatmentArm

# Clear cache of a study object
cvd408$clearCache()

# Connect to the NYVAC durability comparison group
# https://dataspace.cavd.org/cds/CAVD/app.view#group/groupsummary/220
nyvac <- con$getGroup(220)

# Connect to all studies
cvd <- con$getStudy("")

# Refresh the study object to update available datasets and treatment info
cvd$refresh()
}

}
\seealso{
\code{\link{connectDS}} \code{\link{DataSpaceConnection}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{study}}{A character. The study name.}

\item{\code{config}}{A list. Stores configuration of the connection object such
as URL, path and username.}

\item{\code{availableDatasets}}{A data.table. The table of datasets available in
the \code{DataSpaceStudy} object.}

\item{\code{cache}}{A list. Stores the data to avoid downloading the same tables
multiple times.}

\item{\code{dataDir}}{A character. Default directory for storing nonstandard
datasets. Set with \code{setDataDir(dataDir)}.}

\item{\code{treatmentArm}}{A data.table. The table of treatment arm
information for the connected study. Not available for all study
connection.}

\item{\code{group}}{A character. The group name.}

\item{\code{studyInfo}}{A list. Stores the information about the study.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{DataSpaceStudy$new()}}
\item \href{#method-print}{\code{DataSpaceStudy$print()}}
\item \href{#method-getDataset}{\code{DataSpaceStudy$getDataset()}}
\item \href{#method-clearCache}{\code{DataSpaceStudy$clearCache()}}
\item \href{#method-getDatasetDescription}{\code{DataSpaceStudy$getDatasetDescription()}}
\item \href{#method-setDataDir}{\code{DataSpaceStudy$setDataDir()}}
\item \href{#method-refresh}{\code{DataSpaceStudy$refresh()}}
\item \href{#method-clone}{\code{DataSpaceStudy$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize \code{DataSpaceStudy} class.
See \code{\link{DataSpaceConnection}}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceStudy$new(study = NULL, config = NULL, group = NULL, studyInfo = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{study}}{A character. Name of the study to retrieve.}

\item{\code{config}}{A list. Stores configuration of the connection object such
as URL, path and username.}

\item{\code{group}}{An integer. ID of the group to retrieve.}

\item{\code{studyInfo}}{A list. Stores the information about the study.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
Print \code{DataSpaceStudy} class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceStudy$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getDataset"></a>}}
\if{latex}{\out{\hypertarget{method-getDataset}{}}}
\subsection{Method \code{getDataset()}}{
Get a dataset from the connection.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceStudy$getDataset(
  datasetName,
  mergeExtra = FALSE,
  colFilter = NULL,
  reload = FALSE,
  outputDir = NULL,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{datasetName}}{A character. Name of the dataset to retrieve.
Accepts the value in either the "name" or "label" field from \code{availableDatasets}.}

\item{\code{mergeExtra}}{A logical. If set to TRUE, merge extra information.
Ignored for non-integrated datasets.}

\item{\code{colFilter}}{A matrix. A filter as returned by Rlabkey's
\code{\link[Rlabkey]{makeFilter}}.}

\item{\code{reload}}{A logical. If set to TRUE, download the dataset, whether
a cached version exist or not.}

\item{\code{outputDir}}{A character. Optional, specifies directory to download
nonstandard datasets. If \code{NULL}, data will be downloaded to
\code{dataDir}, set with \code{setDataDir(dataDir)}. If \code{dataDir}
is not set, and \code{outputDir} is \code{NULL}, a tmp directory will be
used.}

\item{\code{...}}{Extra arguments to be passed to
\code{\link[Rlabkey]{labkey.selectRows}}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clearCache"></a>}}
\if{latex}{\out{\hypertarget{method-clearCache}{}}}
\subsection{Method \code{clearCache()}}{
Clear \code{cache}. Remove downloaded datasets.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceStudy$clearCache()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getDatasetDescription"></a>}}
\if{latex}{\out{\hypertarget{method-getDatasetDescription}{}}}
\subsection{Method \code{getDatasetDescription()}}{
Get variable information.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceStudy$getDatasetDescription(datasetName, outputDir = NULL)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{datasetName}}{A character. Name of the dataset to retrieve.
Accepts the value in either the "name" or "label" field from \code{availableDatasets}.}

\item{\code{outputDir}}{A character. Directory path.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-setDataDir"></a>}}
\if{latex}{\out{\hypertarget{method-setDataDir}{}}}
\subsection{Method \code{setDataDir()}}{
Set default directory to download non-integrated datasets. If no
\code{dataDir} is set, a tmp directory will be used.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceStudy$setDataDir(dataDir)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dataDir}}{A character. Directory path.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-refresh"></a>}}
\if{latex}{\out{\hypertarget{method-refresh}{}}}
\subsection{Method \code{refresh()}}{
Refresh the study object to update available datasets and treatment info.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceStudy$refresh()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceStudy$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/DataSpaceConnection.R
\name{DataSpaceConnection}
\alias{DataSpaceConnection}
\title{The DataSpaceConnection class}
\description{
The DataSpaceConnection class

The DataSpaceConnection class
}
\section{Constructor}{

\code{\link{connectDS}}
}

\examples{
\dontrun{
# Create a connection (Initiate a DataSpaceConnection object)
con <- connectDS()
con

# Connect to cvd408
# https://dataspace.cavd.org/cds/CAVD/app.view#learn/learn/Study/cvd408?q=408
cvd408 <- con$getStudy("cvd408")

# Connect to all studies
cvd <- con$getStudy("cvd408")

# Connect to the NYVAC durability comparison group
# https://dataspace.cavd.org/cds/CAVD/app.view#group/groupsummary/220
nyvac <- con$getGroup(220)

# Refresh the connection object to update available studies and groups
con$refresh()
}

}
\seealso{
\code{\link{connectDS}} \code{\link{DataSpaceR-package}}
}
\section{Active bindings}{
\if{html}{\out{<div class="r6-active-bindings">}}
\describe{
\item{\code{config}}{A list. Stores configuration of the connection object such as
URL, path and username.}

\item{\code{availableStudies}}{A data.table. The table of available studies.}

\item{\code{availableGroups}}{A data.table. The table of available groups.}

\item{\code{availablePublications}}{A data.table. The table of available
publications.}

\item{\code{mabGridSummary}}{A data.table. The filtered grid with updated
\code{n_} columns and \code{geometric_mean_curve_ic50}.}

\item{\code{mabGrid}}{A data.table. The filtered mAb grid.}

\item{\code{virusMetadata}}{A data.table. Metadata about all viruses in the
DataSpace.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{DataSpaceConnection$new()}}
\item \href{#method-print}{\code{DataSpaceConnection$print()}}
\item \href{#method-getStudy}{\code{DataSpaceConnection$getStudy()}}
\item \href{#method-getGroup}{\code{DataSpaceConnection$getGroup()}}
\item \href{#method-filterMabGrid}{\code{DataSpaceConnection$filterMabGrid()}}
\item \href{#method-resetMabGrid}{\code{DataSpaceConnection$resetMabGrid()}}
\item \href{#method-getMab}{\code{DataSpaceConnection$getMab()}}
\item \href{#method-downloadPublicationData}{\code{DataSpaceConnection$downloadPublicationData()}}
\item \href{#method-refresh}{\code{DataSpaceConnection$refresh()}}
\item \href{#method-clone}{\code{DataSpaceConnection$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Initialize a \code{DataSpaceConnection} object.
See \code{\link{connectDS}}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$new(
  login = NULL,
  password = NULL,
  verbose = FALSE,
  onStaging = FALSE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{login}}{A character. Optional argument. If there is no netrc file a
temporary one can be written by passing login and password of an active
DataSpace account.}

\item{\code{password}}{A character. Optional. The password for the selected
login.}

\item{\code{verbose}}{A logical. Whether to print the extra details for
troubleshooting.}

\item{\code{onStaging}}{A logical. Whether to connect to the staging server
instead of the production server.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new `DataSpaceConnection` object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
Print the \code{DataSpaceConnection} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getStudy"></a>}}
\if{latex}{\out{\hypertarget{method-getStudy}{}}}
\subsection{Method \code{getStudy()}}{
Create a \code{\link{DataSpaceStudy}} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$getStudy(studyName)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{studyName}}{A character. Name of the study to retrieve.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getGroup"></a>}}
\if{latex}{\out{\hypertarget{method-getGroup}{}}}
\subsection{Method \code{getGroup()}}{
Create a \code{\link{DataSpaceStudy}} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$getGroup(groupId)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{groupId}}{An integer. ID of the group to retrieve.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-filterMabGrid"></a>}}
\if{latex}{\out{\hypertarget{method-filterMabGrid}{}}}
\subsection{Method \code{filterMabGrid()}}{
Filter rows in the mAb grid by specifying the values to keep in the
columns found in the \code{mabGrid} field. It takes the column and the
values and filters the underlying tables.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$filterMabGrid(using, value)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{using}}{A character. Name of the column to filter.}

\item{\code{value}}{A character vector. Values to keep in the mAb grid.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-resetMabGrid"></a>}}
\if{latex}{\out{\hypertarget{method-resetMabGrid}{}}}
\subsection{Method \code{resetMabGrid()}}{
Reset the mAb grid to the unfiltered state.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$resetMabGrid()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-getMab"></a>}}
\if{latex}{\out{\hypertarget{method-getMab}{}}}
\subsection{Method \code{getMab()}}{
Create a \code{\link{DataSpaceMab}} object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$getMab()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-downloadPublicationData"></a>}}
\if{latex}{\out{\hypertarget{method-downloadPublicationData}{}}}
\subsection{Method \code{downloadPublicationData()}}{
Download publication data for a chosen publication.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$downloadPublicationData(
  publicationId,
  outputDir = getwd(),
  unzip = TRUE,
  verbose = TRUE
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{publicationId}}{A character/integer. ID for the publication to
download data for.}

\item{\code{outputDir}}{A character. Path to directory to download publication
data.}

\item{\code{unzip}}{A logical. If TRUE, unzip publication data to outputDir.}

\item{\code{verbose}}{A logical. Default TRUE.}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-refresh"></a>}}
\if{latex}{\out{\hypertarget{method-refresh}{}}}
\subsection{Method \code{refresh()}}{
Refresh the connection object to update available studies and groups.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$refresh()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{DataSpaceConnection$clone(deep = FALSE)}\if{html}{\out{</div>}}
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
% Please edit documentation in R/helpers.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{makeFilter}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{Rlabkey}{\code{\link[Rlabkey]{makeFilter}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netrc.R
\name{writeNetrc}
\alias{writeNetrc}
\title{Write a netrc file}
\usage{
writeNetrc(
  login,
  password,
  netrcFile = NULL,
  onStaging = FALSE,
  overwrite = FALSE
)
}
\arguments{
\item{login}{A character. Email address used for logging in on DataSpace.}

\item{password}{A character. Password associated with the login.}

\item{netrcFile}{A character. Credentials will be written into that file.
If left NULL, netrc will be written into a temporary file.}

\item{onStaging}{A logical. Whether to connect to the staging server instead
of the production server.}

\item{overwrite}{A logical. Whether to overwrite the existing netrc file.}
}
\value{
A character vector containing the netrc file path
}
\description{
Write a netrc file that is valid for accessing DataSpace.
}
\details{
The database is accessed with the user's credentials.
A netrc file storing login and password information is required.
See \href{https://docs.ropensci.org/DataSpaceR/}{here}
for instruction on how to register and set DataSpace credential.
By default \code{curl} will look for the file in your home directory.
}
\examples{
# First, create an account in the DataSpace App and read the terms of use
# Next, create a netrc file using writeNetrc()
writeNetrc(
  login = "dataspaceuser@email.com",
  password = "yourSecretPassword"
)
# Specify `netrcFile = getNetrcPath()` to write netrc in the default path
}
\seealso{
\code{\link{connectDS}} \code{\link{checkNetrc}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/netrc.R
\name{getNetrcPath}
\alias{getNetrcPath}
\title{Get a default netrc file path}
\usage{
getNetrcPath()
}
\value{
A character vector containing the default netrc file path
}
\description{
Get a default netrc file path
}
\examples{
getNetrcPath()
}
